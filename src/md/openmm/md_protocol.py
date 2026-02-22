from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path
import time

import numpy as np

from .ligand import build_forcefield, load_ligand_molecule
from .platform import get_platform
from .require import require_module

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class MDProtocolConfig:
    temperature_start_k: float = 10.0
    temperature_target_k: float = 300.0
    pressure_bar: float = 1.0
    timestep_fs: float = 2.0
    minimization_stage1_steps: int = 1000
    minimization_stage2_steps: int = 1000
    minimization_unrestrained_steps: int = 5000
    heating_ps: float = 25.0
    # Canonical panel production length used in this repo's doravirine resistance panel.
    production_ns: float = 10.0
    report_interval_steps: int = 2000
    solvent_padding_nm: float = 1.0
    ionic_strength_molar: float = 0.15
    ca_restraint_k1_kcal_mol_a2: float = 50.0
    ca_restraint_k2_kcal_mol_a2: float = 10.0
    analysis_report_interval_steps: int | None = None


@dataclass(frozen=True)
class MDRunResult:
    total_steps: int
    heating_steps: int
    production_steps: int
    production_steps_completed: int
    elapsed_seconds: float
    resumed_from_checkpoint: bool


def _min_ligand_protein_distance_from_topology_positions(topology, positions, ligand_resname: str) -> float:
    unit = require_module("openmm.unit")
    pos_nm = np.array([p.value_in_unit(unit.nanometer) for p in positions], dtype=float)
    lig_idx: list[int] = []
    prot_idx: list[int] = []
    for atom in topology.atoms():
        if atom.residue.name == ligand_resname:
            lig_idx.append(atom.index)
        elif atom.residue.name in {"HOH", "WAT"}:
            continue
        elif atom.element is not None and atom.element.symbol.upper() in {"NA", "CL", "K"}:
            continue
        else:
            prot_idx.append(atom.index)
    if not lig_idx or not prot_idx:
        return float("nan")
    lig = pos_nm[np.array(lig_idx, dtype=int)]
    prot = pos_nm[np.array(prot_idx, dtype=int)]
    d2 = ((lig[:, None, :] - prot[None, :, :]) ** 2).sum(axis=2)
    return float(np.sqrt(d2.min()) * 10.0)


class _StrippedDCDReporter:
    """DCD reporter that writes only a subset of atoms (protein + ligand).

    Uses openmm.app.DCDFile directly — compatible with OpenMM 8.1.1 which
    lacks the ``atomSubset`` parameter on ``DCDReporter``.
    """

    def __init__(
        self,
        file_path: Path,
        stripped_topology,
        timestep_ps: float,
        interval: int,
        atom_indices: list[int],
        append: bool = False,
        first_step: int = 0,
    ):
        app = require_module("openmm.app")
        unit = require_module("openmm.unit")

        self._interval = int(interval)
        self._atom_indices = atom_indices

        if append and file_path.exists():
            mode = "r+b"
        else:
            mode = "wb"
        self._handle = open(file_path, mode)
        # Always pass the per-frame (not per-integration-step) time step and do
        # NOT pass `interval` to DCDFile.  When `interval` is passed, OpenMM's
        # DCDFile writes nsavc=1 and a garbage DELTA field, which causes
        # downstream readers (MDAnalysis, VMD) to report dt=1.0 ps regardless of
        # the true frame spacing.  Passing dt_frame directly writes the correct
        # DELTA into the DCD header.
        dt_frame = timestep_ps * interval * unit.picoseconds
        try:
            self._dcd = app.DCDFile(
                self._handle,
                stripped_topology,
                dt_frame,
                firstStep=int(first_step),
                append=bool(append),
            )
        except TypeError:
            # Older OpenMM builds lacking firstStep/append kwargs.
            self._dcd = app.DCDFile(self._handle, stripped_topology, dt_frame)

    def describeNextReport(self, simulation):
        steps_done = simulation.currentStep
        steps_left = self._interval - (steps_done % self._interval)
        return (steps_left, True, False, False, False, None)

    def report(self, simulation, state):
        unit = require_module("openmm.unit")
        openmm = require_module("openmm")
        # Convert to plain Vec3 values in nm; this avoids nested Quantity objects
        # that can trigger type errors in OpenMM's DCD writer on some runtimes.
        pos_nm = state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
        subset = [
            openmm.Vec3(float(pos_nm[i, 0]), float(pos_nm[i, 1]), float(pos_nm[i, 2]))
            for i in self._atom_indices
        ]
        self._dcd.writeModel(subset)

    def close(self):
        if hasattr(self, "_handle") and self._handle and not self._handle.closed:
            self._handle.close()

    def __del__(self):
        self.close()


def _write_stripped_topology_pdb(topology, positions, output_path: Path) -> None:
    """Write a PDB containing only solute atoms (no water/ions)."""
    app = require_module("openmm.app")
    solvent_resnames = {"HOH", "WAT"}
    ion_elements = {"NA", "CL", "K"}
    modeller = app.Modeller(topology, positions)
    to_delete = []
    for res in modeller.topology.residues():
        if res.name in solvent_resnames:
            to_delete.append(res)
            continue
        atoms = list(res.atoms())
        if len(atoms) == 1 and atoms[0].element is not None and atoms[0].element.symbol.upper() in ion_elements:
            to_delete.append(res)
    modeller.delete(to_delete)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as handle:
        app.PDBFile.writeFile(modeller.topology, modeller.positions, handle)


def _solute_atom_indices(topology) -> list[int]:
    """Return sorted indices of all non-water, non-ion atoms (protein + ligand + cofactors)."""
    ion_elements = {"NA", "CL", "K"}
    solvent_resnames = {"HOH", "WAT"}
    indices: list[int] = []
    for res in topology.residues():
        if res.name in solvent_resnames:
            continue
        atoms = list(res.atoms())
        if len(atoms) == 1 and atoms[0].element is not None and atoms[0].element.symbol.upper() in ion_elements:
            continue
        for atom in atoms:
            indices.append(atom.index)
    indices.sort()
    return indices


def _add_ca_restraint_force(system, topology, reference_positions, k_kj_mol_nm2: float) -> int:
    openmm = require_module("openmm")
    unit = require_module("openmm.unit")

    force = openmm.CustomExternalForce("k*((x-x0)^2 + (y-y0)^2 + (z-z0)^2)")
    force.addGlobalParameter("k", float(k_kj_mol_nm2))
    force.addPerParticleParameter("x0")
    force.addPerParticleParameter("y0")
    force.addPerParticleParameter("z0")

    for atom in topology.atoms():
        if atom.name == "CA" and atom.residue.name not in {"HOH", "WAT"}:
            pos = reference_positions[atom.index].value_in_unit(unit.nanometer)
            force.addParticle(atom.index, [float(pos[0]), float(pos[1]), float(pos[2])])

    return system.addForce(force)


def prepare_md_assets(
    minimized_pdb_path: Path,
    ligand_resname: str,
    ligand_sdf: Path,
    topology_pdb_path: Path,
    system_xml_path: Path,
    config: MDProtocolConfig | None = None,
) -> None:
    """Prepare solvated explicit-MD assets for Sherlock execution."""
    app = require_module("openmm.app")
    openmm = require_module("openmm")
    unit = require_module("openmm.unit")
    platform, properties = get_platform()

    cfg = config or MDProtocolConfig()

    with open(minimized_pdb_path, "r") as handle:
        pdb = app.PDBFile(handle)

    ligand = load_ligand_molecule(ligand_sdf)
    forcefield = build_forcefield([ligand])

    modeller = app.Modeller(pdb.topology, pdb.positions)
    modeller.addSolvent(
        forcefield,
        model="tip3p",
        padding=cfg.solvent_padding_nm * unit.nanometer,
        ionicStrength=cfg.ionic_strength_molar * unit.molar,
    )

    system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=app.PME,
        nonbondedCutoff=1.0 * unit.nanometer,
        constraints=app.HBonds,
    )

    # Quick local minimization to avoid startup instabilities on cluster.
    integrator = openmm.LangevinMiddleIntegrator(
        cfg.temperature_target_k * unit.kelvin,
        1.0 / unit.picosecond,
        cfg.timestep_fs * unit.femtoseconds,
    )
    simulation = app.Simulation(modeller.topology, system, integrator, platform, properties)
    simulation.context.setPositions(modeller.positions)
    simulation.minimizeEnergy(maxIterations=500)
    pos = simulation.context.getState(getPositions=True).getPositions()
    del simulation  # Release CUDA context before any subsequent simulation.

    min_dist = _min_ligand_protein_distance_from_topology_positions(
        modeller.topology,
        pos,
        ligand_resname,
    )
    if np.isfinite(min_dist) and min_dist > 15.0:
        raise ValueError(
            f"Prepared complex appears unbound (min ligand-protein distance {min_dist:.2f} Å > 15 Å)."
        )

    topology_pdb_path.parent.mkdir(parents=True, exist_ok=True)
    system_xml_path.parent.mkdir(parents=True, exist_ok=True)

    with open(topology_pdb_path, "w") as handle:
        app.PDBFile.writeFile(modeller.topology, pos, handle)
    system_xml_path.write_text(openmm.XmlSerializer.serialize(system))


def prepare_apo_md_assets(
    minimized_pdb_path: Path,
    ligand_resname: str,
    topology_pdb_path: Path,
    system_xml_path: Path,
    config: MDProtocolConfig | None = None,
) -> None:
    """Prepare solvated explicit-MD assets for an apo (ligand-free) system.

    Strips the ligand from the minimized PDB, builds an amber-only forcefield,
    solvates, and writes topology PDB + serialized XML.  The resulting assets
    are drop-in replacements for the holo ones and are executed by the same
    ``run_prepared_md`` / ``src.md.worker`` pathway.
    """
    app = require_module("openmm.app")
    openmm = require_module("openmm")
    unit = require_module("openmm.unit")
    platform, properties = get_platform()

    cfg = config or MDProtocolConfig()

    with open(minimized_pdb_path, "r") as handle:
        pdb = app.PDBFile(handle)

    # Strip ligand from topology/positions before solvation.
    modeller = app.Modeller(pdb.topology, pdb.positions)
    ligand_residues = [res for res in modeller.topology.residues() if res.name == ligand_resname]
    if ligand_residues:
        modeller.delete(ligand_residues)

    forcefield = app.ForceField(
        "amber14/protein.ff14SB.xml",
        "amber14/DNA.bsc1.xml",
        "amber14/tip3p.xml",
    )

    modeller.addSolvent(
        forcefield,
        model="tip3p",
        padding=cfg.solvent_padding_nm * unit.nanometer,
        ionicStrength=cfg.ionic_strength_molar * unit.molar,
    )

    system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=app.PME,
        nonbondedCutoff=1.0 * unit.nanometer,
        constraints=app.HBonds,
    )

    integrator = openmm.LangevinMiddleIntegrator(
        cfg.temperature_target_k * unit.kelvin,
        1.0 / unit.picosecond,
        cfg.timestep_fs * unit.femtoseconds,
    )
    simulation = app.Simulation(modeller.topology, system, integrator, platform, properties)
    simulation.context.setPositions(modeller.positions)
    simulation.minimizeEnergy(maxIterations=500)
    pos = simulation.context.getState(getPositions=True).getPositions()
    del simulation

    topology_pdb_path.parent.mkdir(parents=True, exist_ok=True)
    system_xml_path.parent.mkdir(parents=True, exist_ok=True)

    with open(topology_pdb_path, "w") as handle:
        app.PDBFile.writeFile(modeller.topology, pos, handle)
    system_xml_path.write_text(openmm.XmlSerializer.serialize(system))


def run_prepared_md(
    prepared_topology_pdb: Path,
    prepared_system_xml: Path,
    final_pdb_path: Path,
    state_csv_path: Path | None = None,
    config: MDProtocolConfig | None = None,
    analysis_dcd_path: Path | None = None,
    analysis_topology_pdb_path: Path | None = None,
    checkpoint_path: Path | None = None,
    checkpoint_interval_steps: int = 5000,
    resume_from_checkpoint: bool = True,
) -> MDRunResult:
    """Run minimization -> heating -> NPT production MD from prepared assets."""
    app = require_module("openmm.app")
    openmm = require_module("openmm")
    unit = require_module("openmm.unit")

    cfg = config or MDProtocolConfig()

    with open(prepared_topology_pdb, "r") as handle:
        pdb = app.PDBFile(handle)
    base_system = openmm.XmlSerializer.deserialize(prepared_system_xml.read_text())

    platform, properties = get_platform()
    logger.info("Selected OpenMM platform: %s (properties: %s)", platform.getName(), properties or "default")

    allow_fallback = str(__import__("os").environ.get("OPENMM_ALLOW_FALLBACK", "0")).strip() in {"1", "true", "TRUE", "yes", "YES"}

    def _make_simulation(topology, system, integrator):
        try:
            sim = app.Simulation(topology, system, integrator, platform, properties)
            logger.info("Created simulation on platform: %s", sim.context.getPlatform().getName())
            return sim
        except Exception as exc:
            if not allow_fallback:
                raise RuntimeError(
                    f"Failed to create simulation on {platform.getName()}: {exc}. "
                    f"This may indicate GPU memory exhaustion (e.g. MIG partitions) "
                    f"or no compatible device. Set OPENMM_ALLOW_FALLBACK=1 to try CPU fallback."
                ) from exc
            logger.warning(
                "Failed to initialize platform %s (%s). OPENMM_ALLOW_FALLBACK=1 so falling back to default platform.",
                platform.getName(),
                exc,
            )
            return app.Simulation(topology, system, integrator)
    t0 = time.perf_counter()
    production_steps_target = max(1, int(round((cfg.production_ns * 1_000_000.0) / cfg.timestep_fs)))
    heating_steps = 0
    resumed_from_ckpt = False

    def _create_production_simulation():
        prod_system = openmm.XmlSerializer.deserialize(prepared_system_xml.read_text())
        barostat = openmm.MonteCarloBarostat(
            cfg.pressure_bar * unit.bar,
            cfg.temperature_target_k * unit.kelvin,
            25,
        )
        prod_system.addForce(barostat)
        prod_integrator = openmm.LangevinMiddleIntegrator(
            cfg.temperature_target_k * unit.kelvin,
            1.0 / unit.picosecond,
            cfg.timestep_fs * unit.femtoseconds,
        )
        return _make_simulation(pdb.topology, prod_system, prod_integrator)

    prod = None
    if checkpoint_path is not None and resume_from_checkpoint and checkpoint_path.exists():
        try:
            prod = _create_production_simulation()
            prod.loadCheckpoint(str(checkpoint_path))
            resumed_from_ckpt = True
            logger.info(
                "Resumed production from checkpoint %s at step %d",
                checkpoint_path,
                int(prod.currentStep),
            )
        except Exception as exc:
            logger.warning("Failed to load checkpoint %s (%s). Starting from scratch.", checkpoint_path, exc)
            resumed_from_ckpt = False
            if prod is not None:
                del prod
                prod = None

    if not resumed_from_ckpt:
        # Add C-alpha positional restraints for early minimization stages.
        base_system = openmm.XmlSerializer.deserialize(prepared_system_xml.read_text())
        k_conv = 418.4  # kcal/mol/Å^2 -> kJ/mol/nm^2
        restraint_force_idx = _add_ca_restraint_force(
            base_system,
            pdb.topology,
            pdb.positions,
            cfg.ca_restraint_k1_kcal_mol_a2 * k_conv,
        )
        restraint_force = base_system.getForce(restraint_force_idx)

        integrator = openmm.LangevinMiddleIntegrator(
            cfg.temperature_target_k * unit.kelvin,
            1.0 / unit.picosecond,
            cfg.timestep_fs * unit.femtoseconds,
        )
        sim = _make_simulation(pdb.topology, base_system, integrator)
        sim.context.setPositions(pdb.positions)

        sim.minimizeEnergy(maxIterations=max(1, int(cfg.minimization_stage1_steps)))
        restraint_force.setGlobalParameterDefaultValue(0, cfg.ca_restraint_k2_kcal_mol_a2 * k_conv)
        sim.context.setParameter("k", cfg.ca_restraint_k2_kcal_mol_a2 * k_conv)
        sim.minimizeEnergy(maxIterations=max(1, int(cfg.minimization_stage2_steps)))
        restraint_force.setGlobalParameterDefaultValue(0, 0.0)
        sim.context.setParameter("k", 0.0)
        sim.minimizeEnergy(maxIterations=max(1, int(cfg.minimization_unrestrained_steps)))

        # Heating: 10K -> target temperature (NVT)
        heating_steps = max(1, int(round((cfg.heating_ps * 1000.0) / cfg.timestep_fs)))
        n_chunks = 25
        chunk_steps = max(1, heating_steps // n_chunks)
        for i in range(n_chunks):
            frac = (i + 1) / float(n_chunks)
            temp = cfg.temperature_start_k + frac * (cfg.temperature_target_k - cfg.temperature_start_k)
            integrator.setTemperature(temp * unit.kelvin)
            sim.context.setVelocitiesToTemperature(temp * unit.kelvin)
            sim.step(chunk_steps)

        heated_state = sim.context.getState(getPositions=True, getVelocities=True)
        del sim
        logger.info("Released heating simulation context")
        prod = _create_production_simulation()
        prod.context.setPositions(heated_state.getPositions())
        prod.context.setVelocities(heated_state.getVelocities())

    if checkpoint_path is not None:
        checkpoint_path.parent.mkdir(parents=True, exist_ok=True)
        prod.reporters.append(
            app.CheckpointReporter(
                str(checkpoint_path),
                max(1, int(checkpoint_interval_steps)),
            )
        )

    if state_csv_path is not None:
        state_csv_path.parent.mkdir(parents=True, exist_ok=True)
        prod.reporters.append(
            app.StateDataReporter(
                str(state_csv_path),
                max(1, int(cfg.report_interval_steps)),
                step=True,
                potentialEnergy=True,
                temperature=True,
                volume=True,
                density=True,
                speed=True,
                separator=",",
                append=bool(resumed_from_ckpt and state_csv_path.exists() and state_csv_path.stat().st_size > 0),
            )
        )

    # Optional stripped analysis DCD (protein + ligand only, sparse frames).
    analysis_reporter = None
    if analysis_dcd_path is not None and analysis_topology_pdb_path is not None and cfg.analysis_report_interval_steps is not None:
        solute_idx = _solute_atom_indices(pdb.topology)
        logger.info("Analysis DCD: %d solute atoms, interval=%d steps", len(solute_idx), cfg.analysis_report_interval_steps)
        if not analysis_topology_pdb_path.exists():
            _write_stripped_topology_pdb(pdb.topology, pdb.positions, analysis_topology_pdb_path)
        with open(analysis_topology_pdb_path, "r") as handle:
            analysis_topology = app.PDBFile(handle).topology
        analysis_dcd_path.parent.mkdir(parents=True, exist_ok=True)
        timestep_ps = cfg.timestep_fs / 1000.0
        append_dcd = bool(resumed_from_ckpt and analysis_dcd_path.exists() and analysis_dcd_path.stat().st_size > 0)
        analysis_reporter = _StrippedDCDReporter(
            file_path=analysis_dcd_path,
            stripped_topology=analysis_topology,
            timestep_ps=timestep_ps,
            interval=max(1, int(cfg.analysis_report_interval_steps)),
            atom_indices=solute_idx,
            append=append_dcd,
            first_step=int(prod.currentStep),
        )
        prod.reporters.append(analysis_reporter)

    current_step = int(prod.currentStep)
    remaining_steps = max(0, production_steps_target - current_step)
    if remaining_steps > 0:
        prod.step(remaining_steps)
    else:
        logger.info(
            "Checkpoint already at/above target production steps (%d >= %d); skipping integration.",
            current_step,
            production_steps_target,
        )

    if analysis_reporter is not None:
        analysis_reporter.close()

    final_state = prod.context.getState(getPositions=True)
    final_pdb_path.parent.mkdir(parents=True, exist_ok=True)
    with open(final_pdb_path, "w") as handle:
        app.PDBFile.writeFile(pdb.topology, final_state.getPositions(), handle)

    elapsed = time.perf_counter() - t0
    final_prod_steps = int(prod.currentStep)
    return MDRunResult(
        total_steps=int(heating_steps + final_prod_steps),
        heating_steps=int(heating_steps),
        production_steps=int(production_steps_target),
        production_steps_completed=final_prod_steps,
        elapsed_seconds=float(elapsed),
        resumed_from_checkpoint=bool(resumed_from_ckpt),
    )
