from __future__ import annotations

from dataclasses import dataclass
import logging
import time
from pathlib import Path

import numpy as np

from .platform import get_platform
from .require import require_module


@dataclass(frozen=True)
class AlchemicalConfig:
    temperature_k: float = 300.0
    pressure_bar: float = 1.0
    timestep_fs: float = 2.0
    equilibration_steps: int = 10_000
    production_steps: int = 40_000
    sample_interval: int = 200
    solvent_padding_nm: float = 1.0
    ionic_strength_molar: float = 0.15
    lambda_schedule: tuple[float, ...] = (
        1.0,
        0.95,
        0.9,
        0.8,
        0.7,
        0.6,
        0.5,
        0.4,
        0.3,
        0.2,
        0.1,
        0.05,
        0.0,
    )


@dataclass(frozen=True)
class AlchemicalLegResult:
    delta_g_kj_mol: float
    pair_delta_g_kj_mol: tuple[float, ...]


@dataclass(frozen=True)
class AlchemicalResult:
    complex_leg_kj_mol: float
    solvent_leg_kj_mol: float
    binding_delta_g_kj_mol: float
    lambda_schedule: tuple[float, ...]


def _logmeanexp(values: np.ndarray) -> float:
    vmax = np.max(values)
    return float(vmax + np.log(np.mean(np.exp(values - vmax))))


def _bar_delta_f(w_forward: np.ndarray, w_reverse: np.ndarray) -> float:
    if w_forward.size == 0 or w_reverse.size == 0:
        raise ValueError("BAR requires non-empty forward and reverse work arrays.")

    n_f = float(w_forward.size)
    n_r = float(w_reverse.size)
    c = np.log(n_f / n_r)

    def f(delta_f: float) -> float:
        left = np.sum(1.0 / (1.0 + np.exp(w_forward - delta_f - c)))
        right = np.sum(1.0 / (1.0 + np.exp(w_reverse + delta_f + c)))
        return float(left - right)

    lo = min(float(np.min(w_forward)) - 50.0, -float(np.max(w_reverse)) - 50.0)
    hi = max(float(np.max(w_forward)) + 50.0, -float(np.min(w_reverse)) + 50.0)
    f_lo = f(lo)
    f_hi = f(hi)

    if f_lo * f_hi > 0:
        fwd = -_logmeanexp(-w_forward)
        rev = _logmeanexp(-w_reverse)
        return 0.5 * (fwd + rev)

    for _ in range(200):
        mid = 0.5 * (lo + hi)
        f_mid = f(mid)
        if abs(f_mid) < 1e-7:
            return mid
        if f_lo * f_mid <= 0:
            hi = mid
            f_hi = f_mid
        else:
            lo = mid
            f_lo = f_mid
    return 0.5 * (lo + hi)


def _set_lambda(simulation, value: float) -> None:
    simulation.context.setParameter("lambda_electrostatics", value)
    simulation.context.setParameter("lambda_sterics", value)


def _collect_neighbor_delta_u(simulation, current: float, neighbor: float) -> float:
    unit = require_module("openmm.unit")
    _set_lambda(simulation, current)
    e_current = (
        simulation.context.getState(getEnergy=True)
        .getPotentialEnergy()
        .value_in_unit(unit.kilojoule_per_mole)
    )
    _set_lambda(simulation, neighbor)
    e_neighbor = (
        simulation.context.getState(getEnergy=True)
        .getPotentialEnergy()
        .value_in_unit(unit.kilojoule_per_mole)
    )
    _set_lambda(simulation, current)
    return float(e_neighbor - e_current)


def _run_leg_from_simulation(simulation, config: AlchemicalConfig) -> AlchemicalLegResult:
    unit = require_module("openmm.unit")
    lambdas = tuple(config.lambda_schedule)
    beta = 1.0 / (
        unit.MOLAR_GAS_CONSTANT_R
        * config.temperature_k
        * unit.kelvin
    ).value_in_unit(unit.kilojoule_per_mole)
    forward = [None] * (len(lambdas) - 1)
    reverse = [None] * (len(lambdas) - 1)

    if config.production_steps < config.sample_interval:
        raise ValueError("production_steps must be >= sample_interval.")
    samples = max(1, config.production_steps // config.sample_interval)

    start_wall = None
    for i, lam in enumerate(lambdas):
        if start_wall is None:
            start_wall = time.perf_counter()
        _set_lambda(simulation, lam)
        simulation.step(config.equilibration_steps)

        to_next = []
        to_prev = []
        for _ in range(samples):
            simulation.step(config.sample_interval)
            if i < len(lambdas) - 1:
                delta_u = _collect_neighbor_delta_u(simulation, lam, lambdas[i + 1])
                to_next.append(beta * delta_u)
            if i > 0:
                delta_u = _collect_neighbor_delta_u(simulation, lam, lambdas[i - 1])
                to_prev.append(beta * delta_u)

        if i < len(lambdas) - 1:
            forward[i] = np.array(to_next, dtype=float)
        if i > 0:
            reverse[i - 1] = np.array(to_prev, dtype=float)

        elapsed = time.perf_counter() - start_wall
        done = i + 1
        total = len(lambdas)
        per_window = elapsed / done if done else 0.0
        remaining = per_window * (total - done)
        logging.info(
            "Lambda %d/%d (%.2f) done. Elapsed %.1fs, ETA %.1fs",
            done,
            total,
            lam,
            elapsed,
            remaining,
        )

    pair_delta_f = []
    for i in range(len(lambdas) - 1):
        w_f = forward[i]
        w_r = reverse[i]
        if w_f is None or w_r is None:
            raise RuntimeError(f"Missing BAR work arrays for lambda pair {i}-{i + 1}.")
        pair_delta_f.append(_bar_delta_f(w_f, w_r))

    pair_delta_g = [float(df / beta) for df in pair_delta_f]
    return AlchemicalLegResult(
        delta_g_kj_mol=float(np.sum(pair_delta_g)),
        pair_delta_g_kj_mol=tuple(pair_delta_g),
    )


def _build_alchemical_system(topology, forcefield, ligand_resname: str, config: AlchemicalConfig):
    app = require_module("openmm.app")
    openmm = require_module("openmm")
    unit = require_module("openmm.unit")
    alchemy = require_module("openmmtools.alchemy")

    ligand_atom_indices = [
        atom.index for atom in topology.atoms() if atom.residue.name == ligand_resname
    ]
    if not ligand_atom_indices:
        raise ValueError(f"Ligand {ligand_resname} not found for alchemical leg.")

    base_system = forcefield.createSystem(
        topology,
        nonbondedMethod=app.PME,
        nonbondedCutoff=1.0 * unit.nanometer,
        constraints=app.HBonds,
        rigidWater=True,
    )
    if config.pressure_bar > 0:
        base_system.addForce(
            openmm.MonteCarloBarostat(
                config.pressure_bar * unit.bar,
                config.temperature_k * unit.kelvin,
                25,
            )
        )

    region = alchemy.AlchemicalRegion(
        alchemical_atoms=ligand_atom_indices,
        annihilate_electrostatics=True,
        annihilate_sterics=True,
    )
    factory = alchemy.AbsoluteAlchemicalFactory()
    return factory.create_alchemical_system(base_system, region)


def _run_leg(
    topology,
    positions,
    forcefield,
    ligand_resname: str,
    config: AlchemicalConfig,
) -> AlchemicalLegResult:
    app = require_module("openmm.app")
    openmm = require_module("openmm")
    unit = require_module("openmm.unit")
    platform, properties = get_platform()
    alchemical_system = _build_alchemical_system(
        topology=topology,
        forcefield=forcefield,
        ligand_resname=ligand_resname,
        config=config,
    )

    integrator = openmm.LangevinMiddleIntegrator(
        config.temperature_k * unit.kelvin,
        1.0 / unit.picosecond,
        config.timestep_fs * unit.femtoseconds,
    )
    simulation = app.Simulation(
        topology, alchemical_system, integrator, platform, properties
    )
    simulation.context.setPositions(positions)
    simulation.context.setVelocitiesToTemperature(config.temperature_k * unit.kelvin)
    return _run_leg_from_simulation(simulation, config)


def _extract_ligand_only(topology, positions, ligand_resname: str):
    app = require_module("openmm.app")
    modeller = app.Modeller(topology, positions)
    to_delete = [res for res in modeller.topology.residues() if res.name != ligand_resname]
    modeller.delete(to_delete)
    return modeller.topology, modeller.positions


def _solvate(topology, positions, forcefield, config: AlchemicalConfig):
    app = require_module("openmm.app")
    unit = require_module("openmm.unit")
    modeller = app.Modeller(topology, positions)
    modeller.addSolvent(
        forcefield,
        model="tip3p",
        padding=config.solvent_padding_nm * unit.nanometer,
        ionicStrength=config.ionic_strength_molar * unit.molar,
    )
    return modeller.topology, modeller.positions


def compute_alchemical_binding_free_energy(
    topology,
    positions,
    forcefield,
    ligand_resname: str,
    config: AlchemicalConfig | None = None,
) -> AlchemicalResult:
    cfg = config or AlchemicalConfig()

    complex_top, complex_pos = _solvate(topology, positions, forcefield, cfg)
    complex_leg = _run_leg(
        topology=complex_top,
        positions=complex_pos,
        forcefield=forcefield,
        ligand_resname=ligand_resname,
        config=cfg,
    )

    ligand_top, ligand_pos = _extract_ligand_only(topology, positions, ligand_resname)
    solv_lig_top, solv_lig_pos = _solvate(ligand_top, ligand_pos, forcefield, cfg)
    solvent_leg = _run_leg(
        topology=solv_lig_top,
        positions=solv_lig_pos,
        forcefield=forcefield,
        ligand_resname=ligand_resname,
        config=cfg,
    )

    binding = complex_leg.delta_g_kj_mol - solvent_leg.delta_g_kj_mol
    return AlchemicalResult(
        complex_leg_kj_mol=complex_leg.delta_g_kj_mol,
        solvent_leg_kj_mol=solvent_leg.delta_g_kj_mol,
        binding_delta_g_kj_mol=binding,
        lambda_schedule=cfg.lambda_schedule,
    )


def write_alchemical_result(path: Path, result: AlchemicalResult, metadata: dict) -> None:
    import json

    payload = {
        **metadata,
        "complex_leg_kj_mol": result.complex_leg_kj_mol,
        "solvent_leg_kj_mol": result.solvent_leg_kj_mol,
        "binding_delta_g_kj_mol": result.binding_delta_g_kj_mol,
        "lambda_schedule": list(result.lambda_schedule),
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2))


@dataclass(frozen=True)
class SingleLegResult:
    """Result from running a single FEP leg (complex or solvent)."""

    leg: str
    delta_g_kj_mol: float
    pair_delta_g_kj_mol: tuple[float, ...]
    lambda_schedule: tuple[float, ...]


def write_single_leg_result(
    path: Path, result: SingleLegResult, metadata: dict
) -> None:
    """Write a single leg FEP result to JSON."""
    import json

    payload = {
        **metadata,
        "leg": result.leg,
        "delta_g_kj_mol": result.delta_g_kj_mol,
        "pair_delta_g_kj_mol": list(result.pair_delta_g_kj_mol),
        "lambda_schedule": list(result.lambda_schedule),
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2))


def run_single_leg(
    minimized_pdb_path: Path,
    ligand_resname: str,
    ligand_sdf: Path,
    leg: str,
    config: AlchemicalConfig | None = None,
    output_json: Path | None = None,
    metadata: dict | None = None,
) -> SingleLegResult:
    """Run a single FEP leg (complex or solvent) for cluster execution.

    Args:
        minimized_pdb_path: Path to pre-minimized PDB structure.
        ligand_resname: Residue name of the ligand.
        ligand_sdf: Path to ligand SDF file.
        leg: Either "complex" or "solvent".
        config: Alchemical configuration parameters.
        output_json: Optional path to save result as JSON.
        metadata: Optional metadata to include in JSON output.

    Returns:
        SingleLegResult containing the free energy for this leg.
    """
    from .pipeline import build_forcefield, load_ligand_molecule

    if leg not in ("complex", "solvent"):
        raise ValueError(f"leg must be 'complex' or 'solvent', got {leg!r}")

    app = require_module("openmm.app")
    cfg = config or AlchemicalConfig()

    with open(minimized_pdb_path, "r") as handle:
        pdb = app.PDBFile(handle)

    ligand = load_ligand_molecule(ligand_sdf)
    forcefield = build_forcefield([ligand])

    if leg == "complex":
        top, pos = _solvate(pdb.topology, pdb.positions, forcefield, cfg)
    else:
        lig_top, lig_pos = _extract_ligand_only(
            pdb.topology, pdb.positions, ligand_resname
        )
        top, pos = _solvate(lig_top, lig_pos, forcefield, cfg)

    leg_result = _run_leg(
        topology=top,
        positions=pos,
        forcefield=forcefield,
        ligand_resname=ligand_resname,
        config=cfg,
    )

    result = SingleLegResult(
        leg=leg,
        delta_g_kj_mol=leg_result.delta_g_kj_mol,
        pair_delta_g_kj_mol=leg_result.pair_delta_g_kj_mol,
        lambda_schedule=cfg.lambda_schedule,
    )

    if output_json is not None:
        write_single_leg_result(output_json, result, metadata or {})

    return result


def prepare_single_leg_assets(
    minimized_pdb_path: Path,
    ligand_resname: str,
    ligand_sdf: Path,
    leg: str,
    topology_pdb_path: Path,
    system_xml_path: Path,
    config: AlchemicalConfig | None = None,
) -> None:
    """Prepare serialized alchemical assets locally for OpenMM-only cluster execution."""
    from .pipeline import build_forcefield, load_ligand_molecule

    if leg not in ("complex", "solvent"):
        raise ValueError(f"leg must be 'complex' or 'solvent', got {leg!r}")

    app = require_module("openmm.app")
    openmm = require_module("openmm")
    cfg = config or AlchemicalConfig()

    with open(minimized_pdb_path, "r") as handle:
        pdb = app.PDBFile(handle)

    ligand = load_ligand_molecule(ligand_sdf)
    forcefield = build_forcefield([ligand])

    if leg == "complex":
        top, pos = _solvate(pdb.topology, pdb.positions, forcefield, cfg)
    else:
        lig_top, lig_pos = _extract_ligand_only(
            pdb.topology, pdb.positions, ligand_resname
        )
        top, pos = _solvate(lig_top, lig_pos, forcefield, cfg)

    system = _build_alchemical_system(
        topology=top,
        forcefield=forcefield,
        ligand_resname=ligand_resname,
        config=cfg,
    )

    # Minimize solvated alchemical system locally to reduce NaNs on cluster.
    openmm = require_module("openmm")
    unit = require_module("openmm.unit")
    integrator = openmm.LangevinMiddleIntegrator(
        cfg.temperature_k * unit.kelvin,
        1.0 / unit.picosecond,
        cfg.timestep_fs * unit.femtoseconds,
    )
    simulation = app.Simulation(top, system, integrator)
    simulation.context.setPositions(pos)
    simulation.minimizeEnergy()
    pos = simulation.context.getState(getPositions=True).getPositions()

    topology_pdb_path.parent.mkdir(parents=True, exist_ok=True)
    system_xml_path.parent.mkdir(parents=True, exist_ok=True)
    with open(topology_pdb_path, "w") as handle:
        app.PDBFile.writeFile(top, pos, handle)
    system_xml_path.write_text(openmm.XmlSerializer.serialize(system))


def run_single_leg_prepared(
    prepared_topology_pdb: Path,
    prepared_system_xml: Path,
    leg: str,
    config: AlchemicalConfig | None = None,
    output_json: Path | None = None,
    metadata: dict | None = None,
) -> SingleLegResult:
    """Run a single FEP leg from prebuilt serialized assets (OpenMM-only runtime)."""
    if leg not in ("complex", "solvent"):
        raise ValueError(f"leg must be 'complex' or 'solvent', got {leg!r}")

    app = require_module("openmm.app")
    openmm = require_module("openmm")
    unit = require_module("openmm.unit")
    cfg = config or AlchemicalConfig()
    platform, properties = get_platform()

    with open(prepared_topology_pdb, "r") as handle:
        pdb = app.PDBFile(handle)
    system = openmm.XmlSerializer.deserialize(prepared_system_xml.read_text())
    integrator = openmm.LangevinMiddleIntegrator(
        cfg.temperature_k * unit.kelvin,
        1.0 / unit.picosecond,
        cfg.timestep_fs * unit.femtoseconds,
    )
    simulation = app.Simulation(
        pdb.topology, system, integrator, platform, properties
    )
    simulation.context.setPositions(pdb.positions)
    # Minimize solvated system to avoid NaNs during equilibration.
    simulation.minimizeEnergy()
    simulation.context.setVelocitiesToTemperature(cfg.temperature_k * unit.kelvin)
    # Gentle warmup with a smaller timestep to reduce initial instabilities.
    try:
        original_step = integrator.getStepSize()
        warmup_step = min(original_step, 0.5 * unit.femtoseconds)
        if warmup_step < original_step:
            integrator.setStepSize(warmup_step)
        simulation.step(200)
        if warmup_step < original_step:
            integrator.setStepSize(original_step)
    except Exception:
        simulation.step(100)

    leg_result = _run_leg_from_simulation(simulation, cfg)
    result = SingleLegResult(
        leg=leg,
        delta_g_kj_mol=leg_result.delta_g_kj_mol,
        pair_delta_g_kj_mol=leg_result.pair_delta_g_kj_mol,
        lambda_schedule=cfg.lambda_schedule,
    )

    if output_json is not None:
        write_single_leg_result(output_json, result, metadata or {})

    return result
