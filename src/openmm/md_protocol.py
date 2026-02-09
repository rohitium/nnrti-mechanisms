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
    production_ns: float = 2.0
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
    elapsed_seconds: float


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

    def __init__(self, file_path: Path, topology, timestep_ps: float, interval: int, atom_indices: list[int]):
        app = require_module("openmm.app")
        openmm = require_module("openmm")
        unit = require_module("openmm.unit")

        self._interval = int(interval)
        self._atom_indices = atom_indices

        # Build stripped topology by deleting solvent/ion residues.
        solvent_resnames = {"HOH", "WAT"}
        ion_elements = {"NA", "CL", "K"}
        n_atoms = sum(1 for _ in topology.atoms())
        dummy_pos = [openmm.Vec3(0, 0, 0) * unit.nanometer] * n_atoms
        modeller = app.Modeller(topology, dummy_pos)
        to_delete = []
        for res in modeller.topology.residues():
            if res.name in solvent_resnames:
                to_delete.append(res)
                continue
            atoms = list(res.atoms())
            if len(atoms) == 1 and atoms[0].element is not None and atoms[0].element.symbol.upper() in ion_elements:
                to_delete.append(res)
        modeller.delete(to_delete)

        self._handle = open(file_path, "wb")
        dt = timestep_ps * interval * unit.picoseconds
        self._dcd = app.DCDFile(self._handle, modeller.topology, dt)

    def describeNextReport(self, simulation):
        steps_done = simulation.currentStep
        steps_left = self._interval - (steps_done % self._interval)
        return (steps_left, True, False, False, False, None)

    def report(self, simulation, state):
        unit = require_module("openmm.unit")
        positions = state.getPositions()
        subset = [positions[i] for i in self._atom_indices]
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
    for atom in topology.atoms():
        if atom.residue.name in solvent_resnames:
            continue
        if atom.element is not None and atom.element.symbol.upper() in ion_elements:
            continue
        indices.append(atom.index)
    indices.sort()
    # Solute atoms must be contiguous starting at 0 (addSolvent appends after solute).
    if indices and (indices[0] != 0 or indices[-1] != len(indices) - 1):
        raise RuntimeError(
            f"Solute atom indices are not contiguous [0..{len(indices)-1}]: "
            f"first={indices[0]}, last={indices[-1]}. "
            "This breaks the assumption that addSolvent appends after solute."
        )
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
    simulation = app.Simulation(modeller.topology, system, integrator)
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


def run_prepared_md(
    prepared_topology_pdb: Path,
    prepared_system_xml: Path,
    final_pdb_path: Path,
    state_csv_path: Path | None = None,
    config: MDProtocolConfig | None = None,
    analysis_dcd_path: Path | None = None,
    analysis_topology_pdb_path: Path | None = None,
) -> MDRunResult:
    """Run minimization -> heating -> NPT production MD from prepared assets."""
    app = require_module("openmm.app")
    openmm = require_module("openmm")
    unit = require_module("openmm.unit")

    cfg = config or MDProtocolConfig()

    with open(prepared_topology_pdb, "r") as handle:
        pdb = app.PDBFile(handle)
    base_system = openmm.XmlSerializer.deserialize(prepared_system_xml.read_text())

    # Add C-alpha positional restraints for the early minimization stages.
    k_conv = 418.4  # kcal/mol/Å^2 -> kJ/mol/nm^2
    restraint_force_idx = _add_ca_restraint_force(
        base_system,
        pdb.topology,
        pdb.positions,
        cfg.ca_restraint_k1_kcal_mol_a2 * k_conv,
    )
    restraint_force = base_system.getForce(restraint_force_idx)

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
    integrator = openmm.LangevinMiddleIntegrator(
        cfg.temperature_target_k * unit.kelvin,
        1.0 / unit.picosecond,
        cfg.timestep_fs * unit.femtoseconds,
    )
    sim = _make_simulation(pdb.topology, base_system, integrator)
    sim.context.setPositions(pdb.positions)

    t0 = time.perf_counter()

    sim.minimizeEnergy(maxIterations=max(1, int(cfg.minimization_stage1_steps)))
    restraint_force.setGlobalParameterDefaultValue(0, cfg.ca_restraint_k2_kcal_mol_a2 * k_conv)
    sim.context.setParameter("k", cfg.ca_restraint_k2_kcal_mol_a2 * k_conv)
    sim.minimizeEnergy(maxIterations=max(1, int(cfg.minimization_stage2_steps)))
    restraint_force.setGlobalParameterDefaultValue(0, 0.0)
    sim.context.setParameter("k", 0.0)
    sim.minimizeEnergy(maxIterations=max(1, int(cfg.minimization_unrestrained_steps)))

    # Heating: 10K -> 300K over 25 ps (NVT)
    heating_steps = max(1, int(round((cfg.heating_ps * 1000.0) / cfg.timestep_fs)))
    n_chunks = 25
    chunk_steps = max(1, heating_steps // n_chunks)
    for i in range(n_chunks):
        frac = (i + 1) / float(n_chunks)
        temp = cfg.temperature_start_k + frac * (cfg.temperature_target_k - cfg.temperature_start_k)
        integrator.setTemperature(temp * unit.kelvin)
        sim.context.setVelocitiesToTemperature(temp * unit.kelvin)
        sim.step(chunk_steps)

    state = sim.context.getState(getPositions=True, getVelocities=True)

    # Release the heating CUDA context before creating the production simulation.
    # On memory-limited GPUs (MIG partitions, etc.) two contexts cannot coexist.
    del sim
    logger.info("Released heating simulation context")

    # NPT production simulation.
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
    prod = _make_simulation(pdb.topology, prod_system, prod_integrator)
    prod.context.setPositions(state.getPositions())
    prod.context.setVelocities(state.getVelocities())

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
            )
        )

    # Optional stripped analysis DCD (protein + ligand only, sparse frames).
    analysis_reporter = None
    if analysis_dcd_path is not None and analysis_topology_pdb_path is not None and cfg.analysis_report_interval_steps is not None:
        solute_idx = _solute_atom_indices(pdb.topology)
        logger.info("Analysis DCD: %d solute atoms, interval=%d steps", len(solute_idx), cfg.analysis_report_interval_steps)
        _write_stripped_topology_pdb(pdb.topology, pdb.positions, analysis_topology_pdb_path)
        analysis_dcd_path.parent.mkdir(parents=True, exist_ok=True)
        timestep_ps = cfg.timestep_fs / 1000.0
        analysis_reporter = _StrippedDCDReporter(
            file_path=analysis_dcd_path,
            topology=pdb.topology,
            timestep_ps=timestep_ps,
            interval=max(1, int(cfg.analysis_report_interval_steps)),
            atom_indices=solute_idx,
        )
        prod.reporters.append(analysis_reporter)

    production_steps = max(1, int(round((cfg.production_ns * 1_000_000.0) / cfg.timestep_fs)))
    prod.step(production_steps)

    if analysis_reporter is not None:
        analysis_reporter.close()

    final_state = prod.context.getState(getPositions=True)
    final_pdb_path.parent.mkdir(parents=True, exist_ok=True)
    with open(final_pdb_path, "w") as handle:
        app.PDBFile.writeFile(pdb.topology, final_state.getPositions(), handle)

    elapsed = time.perf_counter() - t0
    return MDRunResult(
        total_steps=int(heating_steps + production_steps),
        heating_steps=int(heating_steps),
        production_steps=int(production_steps),
        elapsed_seconds=float(elapsed),
    )
