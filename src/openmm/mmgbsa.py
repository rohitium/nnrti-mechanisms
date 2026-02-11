from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np

from .ligand import build_forcefield, load_ligand_molecule
from .require import require_module


@dataclass(frozen=True)
class MMGBSAResult:
    n_snapshots: int
    snapshot_indices: tuple[int, ...]
    delta_e_vdw_mean: float
    delta_e_vdw_std: float
    delta_e_vdw_sem: float
    delta_e_elec_mean: float
    delta_e_elec_std: float
    delta_e_elec_sem: float
    delta_g_gb_mean: float
    delta_g_gb_std: float
    delta_g_gb_sem: float
    delta_g_sa_mean: float
    delta_g_sa_std: float
    delta_g_sa_sem: float
    binding_dg_mean: float
    binding_dg_std: float
    binding_dg_sem: float


def _radius_from_symbol(symbol: str) -> float:
    table = {
        "H": 0.12,
        "C": 0.17,
        "N": 0.155,
        "O": 0.152,
        "F": 0.147,
        "P": 0.18,
        "S": 0.18,
        "CL": 0.175,
        "BR": 0.185,
        "I": 0.198,
    }
    return table.get(symbol.upper(), 0.17)


def _build_component_system(topology, forcefield, mode: str):
    app = require_module("openmm.app")
    openmm = require_module("openmm")

    system = forcefield.createSystem(
        topology,
        nonbondedMethod=app.NoCutoff,
        constraints=None,
    )

    for force in system.getForces():
        force.setForceGroup(0)

    nb = None
    for force in system.getForces():
        if isinstance(force, openmm.NonbondedForce):
            nb = force
            break
    if nb is None:
        raise RuntimeError("NonbondedForce not found.")

    if mode == "elec":
        for i in range(nb.getNumParticles()):
            q, sigma, _eps = nb.getParticleParameters(i)
            nb.setParticleParameters(i, q, sigma, 0.0)
        for i in range(nb.getNumExceptions()):
            p1, p2, qprod, sigma, _eps = nb.getExceptionParameters(i)
            nb.setExceptionParameters(i, p1, p2, qprod, sigma, 0.0)
    elif mode == "vdw":
        for i in range(nb.getNumParticles()):
            _q, sigma, eps = nb.getParticleParameters(i)
            nb.setParticleParameters(i, 0.0, sigma, eps)
        for i in range(nb.getNumExceptions()):
            p1, p2, _qprod, sigma, eps = nb.getExceptionParameters(i)
            nb.setExceptionParameters(i, p1, p2, 0.0, sigma, eps)
    elif mode != "full":
        raise ValueError(f"Unknown mode: {mode}")

    nb.setForceGroup(1)
    return system


def _build_gb_system(topology, forcefield, include_sa: bool):
    app = require_module("openmm.app")
    openmm = require_module("openmm")

    system = forcefield.createSystem(
        topology,
        nonbondedMethod=app.NoCutoff,
        constraints=None,
    )
    for force in system.getForces():
        force.setForceGroup(0)

    nb = None
    for force in system.getForces():
        if isinstance(force, openmm.NonbondedForce):
            nb = force
            break
    if nb is None:
        raise RuntimeError("NonbondedForce not found for GB system.")

    gb = openmm.GBSAOBCForce()
    gb.setSolventDielectric(80.0)
    gb.setSoluteDielectric(2.0)
    if include_sa:
        gb.setSurfaceAreaEnergy(2.25936)  # kJ/mol/nm^2
    else:
        gb.setSurfaceAreaEnergy(0.0)

    for atom in topology.atoms():
        q, _sigma, _eps = nb.getParticleParameters(atom.index)
        symbol = atom.element.symbol if atom.element is not None else atom.name[0]
        gb.addParticle(q, _radius_from_symbol(symbol), 1.0)

    gb.setForceGroup(2)
    system.addForce(gb)
    return system


def _make_context(system):
    openmm = require_module("openmm")
    unit = require_module("openmm.unit")
    integrator = openmm.VerletIntegrator(0.001 * unit.picoseconds)
    return openmm.Context(system, integrator)


def _energy_of(context, positions_nm: np.ndarray, force_group: int) -> float:
    unit = require_module("openmm.unit")
    openmm = require_module("openmm")

    pos = [openmm.Vec3(float(x), float(y), float(z)) * unit.nanometer for x, y, z in positions_nm]
    context.setPositions(pos)
    state = context.getState(getEnergy=True, groups=1 << int(force_group))
    return float(state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole))


def _subset_positions(all_positions_nm: np.ndarray, indices: np.ndarray) -> np.ndarray:
    return all_positions_nm[indices]


def _select_snapshot_indices(
    n_frames: int,
    discard_fraction: float,
    n_snapshots: int,
    dt_ps: float | None = None,
    sample_window_ns: float | None = 1.0,
) -> np.ndarray:
    if (
        sample_window_ns is not None
        and sample_window_ns > 0.0
        and dt_ps is not None
        and np.isfinite(dt_ps)
        and dt_ps > 0.0
        and n_frames > 1
    ):
        total_time_ps = float(n_frames - 1) * float(dt_ps)
        window_ps = float(sample_window_ns) * 1000.0
        start_time_ps = max(0.0, total_time_ps - window_ps)
        start = int(np.ceil(start_time_ps / float(dt_ps)))
        start = min(max(0, start), max(0, n_frames - 1))
    else:
        # Backward-compatible fallback when trajectory timing metadata is unavailable.
        start = int(np.floor(max(0.0, min(0.95, discard_fraction)) * n_frames))
        start = min(start, max(0, n_frames - 1))
    available = n_frames - start
    n_take = min(max(1, int(n_snapshots)), max(1, available))
    return np.linspace(start, n_frames - 1, num=n_take, dtype=int)


def _extract_receptor_indices(topology, ligand_resname: str) -> np.ndarray:
    idx: list[int] = []
    for atom in topology.atoms():
        if atom.residue.name != ligand_resname:
            idx.append(atom.index)
    return np.asarray(idx, dtype=int)


def _extract_ligand_indices(topology, ligand_resname: str) -> np.ndarray:
    idx: list[int] = []
    for atom in topology.atoms():
        if atom.residue.name == ligand_resname:
            idx.append(atom.index)
    return np.asarray(idx, dtype=int)


def _make_subtopology(topology, positions, keep_residue_name: str | None, invert: bool):
    app = require_module("openmm.app")
    modeller = app.Modeller(topology, positions)
    to_delete = []
    for res in modeller.topology.residues():
        is_lig = res.name == keep_residue_name
        if invert:
            if is_lig:
                to_delete.append(res)
        else:
            if not is_lig:
                to_delete.append(res)
    modeller.delete(to_delete)
    return modeller.topology


def compute_mmgbsa_from_trajectory(
    minimized_pdb_path: Path,
    trajectory_dcd_path: Path,
    ligand_resname: str,
    ligand_sdf: Path,
    n_snapshots: int = 100,
    discard_fraction: float = 0.25,
    sample_window_ns: float | None = 1.0,
    analysis_topology_pdb_path: Path | None = None,
) -> MMGBSAResult:
    """Compute MM/GBSA-style decomposition from explicit-MD snapshots."""
    import MDAnalysis as mda

    app = require_module("openmm.app")

    with open(minimized_pdb_path, "r") as handle:
        solute_pdb = app.PDBFile(handle)

    ligand = load_ligand_molecule(ligand_sdf)
    forcefield = build_forcefield([ligand])

    complex_top = solute_pdb.topology
    complex_n = sum(1 for _ in complex_top.atoms())

    receptor_idx = _extract_receptor_indices(complex_top, ligand_resname)
    ligand_idx = _extract_ligand_indices(complex_top, ligand_resname)
    if receptor_idx.size == 0 or ligand_idx.size == 0:
        raise ValueError("Could not identify receptor/ligand atoms for MM/GBSA decomposition.")

    receptor_top = _make_subtopology(complex_top, solute_pdb.positions, ligand_resname, invert=True)
    ligand_top = _make_subtopology(complex_top, solute_pdb.positions, ligand_resname, invert=False)

    systems = {
        "complex_vdw": _build_component_system(complex_top, forcefield, "vdw"),
        "complex_elec": _build_component_system(complex_top, forcefield, "elec"),
        "complex_gb_total": _build_gb_system(complex_top, forcefield, include_sa=True),
        "complex_gb_polar": _build_gb_system(complex_top, forcefield, include_sa=False),
        "receptor_vdw": _build_component_system(receptor_top, forcefield, "vdw"),
        "receptor_elec": _build_component_system(receptor_top, forcefield, "elec"),
        "receptor_gb_total": _build_gb_system(receptor_top, forcefield, include_sa=True),
        "receptor_gb_polar": _build_gb_system(receptor_top, forcefield, include_sa=False),
        "ligand_vdw": _build_component_system(ligand_top, forcefield, "vdw"),
        "ligand_elec": _build_component_system(ligand_top, forcefield, "elec"),
        "ligand_gb_total": _build_gb_system(ligand_top, forcefield, include_sa=True),
        "ligand_gb_polar": _build_gb_system(ligand_top, forcefield, include_sa=False),
    }
    contexts = {k: _make_context(v) for k, v in systems.items()}

    # Load trajectory — use stripped topology PDB if provided (solute-only DCD),
    # otherwise fall back to minimized PDB as topology (assumes DCD matches).
    traj_topo = analysis_topology_pdb_path if analysis_topology_pdb_path is not None else minimized_pdb_path
    u = mda.Universe(str(traj_topo), str(trajectory_dcd_path))
    n_frames = len(u.trajectory)
    if n_frames < 1:
        raise ValueError("Empty trajectory for MM/GBSA evaluation.")

    snap_idx = _select_snapshot_indices(
        n_frames=n_frames,
        discard_fraction=discard_fraction,
        n_snapshots=n_snapshots,
        dt_ps=getattr(u.trajectory, "dt", None),
        sample_window_ns=sample_window_ns,
    )

    d_vdw: list[float] = []
    d_elec: list[float] = []
    d_gb: list[float] = []
    d_sa: list[float] = []
    d_tot: list[float] = []

    for frame in snap_idx:
        u.trajectory[int(frame)]
        pos_a = np.asarray(u.atoms.positions, dtype=float)
        if pos_a.shape[0] != complex_n:
            raise ValueError(
                f"Trajectory has {pos_a.shape[0]} atoms but expected {complex_n} (complex_n from minimized PDB)."
            )
        solute_nm = pos_a / 10.0
        rec_nm = _subset_positions(solute_nm, receptor_idx)
        lig_nm = _subset_positions(solute_nm, ligand_idx)

        e_vdw_c = _energy_of(contexts["complex_vdw"], solute_nm, 1)
        e_vdw_r = _energy_of(contexts["receptor_vdw"], rec_nm, 1)
        e_vdw_l = _energy_of(contexts["ligand_vdw"], lig_nm, 1)

        e_elec_c = _energy_of(contexts["complex_elec"], solute_nm, 1)
        e_elec_r = _energy_of(contexts["receptor_elec"], rec_nm, 1)
        e_elec_l = _energy_of(contexts["ligand_elec"], lig_nm, 1)

        gb_t_c = _energy_of(contexts["complex_gb_total"], solute_nm, 2)
        gb_t_r = _energy_of(contexts["receptor_gb_total"], rec_nm, 2)
        gb_t_l = _energy_of(contexts["ligand_gb_total"], lig_nm, 2)

        gb_p_c = _energy_of(contexts["complex_gb_polar"], solute_nm, 2)
        gb_p_r = _energy_of(contexts["receptor_gb_polar"], rec_nm, 2)
        gb_p_l = _energy_of(contexts["ligand_gb_polar"], lig_nm, 2)

        dv = e_vdw_c - e_vdw_r - e_vdw_l
        de = e_elec_c - e_elec_r - e_elec_l
        dgb = gb_p_c - gb_p_r - gb_p_l
        dsa = (gb_t_c - gb_p_c) - (gb_t_r - gb_p_r) - (gb_t_l - gb_p_l)
        dtotal = dv + de + dgb + dsa

        d_vdw.append(float(dv))
        d_elec.append(float(de))
        d_gb.append(float(dgb))
        d_sa.append(float(dsa))
        d_tot.append(float(dtotal))

    def _stats(values: Iterable[float]) -> tuple[float, float, float]:
        arr = np.asarray(list(values), dtype=float)
        mean = float(np.nanmean(arr))
        std = float(np.nanstd(arr, ddof=1)) if arr.size > 1 else 0.0
        sem = float(std / np.sqrt(arr.size)) if arr.size > 0 else float("nan")
        return mean, std, sem

    v_m, v_s, v_se = _stats(d_vdw)
    e_m, e_s, e_se = _stats(d_elec)
    g_m, g_s, g_se = _stats(d_gb)
    s_m, s_s, s_se = _stats(d_sa)
    t_m, t_s, t_se = _stats(d_tot)

    return MMGBSAResult(
        n_snapshots=len(d_tot),
        snapshot_indices=tuple(int(x) for x in snap_idx.tolist()),
        delta_e_vdw_mean=v_m,
        delta_e_vdw_std=v_s,
        delta_e_vdw_sem=v_se,
        delta_e_elec_mean=e_m,
        delta_e_elec_std=e_s,
        delta_e_elec_sem=e_se,
        delta_g_gb_mean=g_m,
        delta_g_gb_std=g_s,
        delta_g_gb_sem=g_se,
        delta_g_sa_mean=s_m,
        delta_g_sa_std=s_s,
        delta_g_sa_sem=s_se,
        binding_dg_mean=t_m,
        binding_dg_std=t_s,
        binding_dg_sem=t_se,
    )
