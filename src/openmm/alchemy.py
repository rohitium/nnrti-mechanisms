from __future__ import annotations

from dataclasses import dataclass
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
    alchemy = require_module("openmmtools.alchemy")
    platform, properties = get_platform()

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
    alchemical_system = factory.create_alchemical_system(base_system, region)

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

    for i, lam in enumerate(lambdas):
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
