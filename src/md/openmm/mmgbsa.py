from __future__ import annotations

from dataclasses import dataclass
import os
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
    #: Absolute (unsubtracted) energies keyed "<phase>_<term>" for phase in
    #: complex/receptor/ligand and term in vdw/elec/gb/sa/total, each mapping to
    #: (mean, std, sem) over the sampled snapshots. The binding terms above are
    #: complex - receptor - ligand of these; kept separately so the raw phase
    #: energies can be inspected without recomputing.
    absolute_terms: tuple[tuple[str, tuple[float, float, float]], ...] = ()


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


def _obc_scale_factors() -> dict[str, float]:
    """Element -> OBC screening factor, from OpenMM's own reference table.

    ``GBSAOBCForce.addParticle`` takes (charge, radius, scaleFactor). The scale
    factor is element-specific in the OBC model (C 0.72, H 0.85, N 0.79, ...);
    passing a constant distorts every Born radius in a burial-dependent way,
    which is exactly what the GB term measures. Elements absent from OpenMM's
    table (e.g. Cl) fall back to the same default OpenMM itself uses.
    """
    from openmm.app.internal.customgbforces import _SCREEN_PARAMETERS

    table: dict[str, float] = {}
    default = _OBC_DEFAULT_SCALE
    for element, values in _SCREEN_PARAMETERS.items():
        if element is None:
            default = float(values[0])
            continue
        table[element.symbol.upper()] = float(values[0])
    table.setdefault("_default", default)
    return table


#: Fallback screening factor for elements OpenMM does not tabulate.
_OBC_DEFAULT_SCALE = 0.8

#: Interior dielectric for the GB term. The MM electrostatics are evaluated in
#: vacuum (NoCutoff NonbondedForce, eps = 1), so the GB term describes a
#: vacuum -> water transfer and eps_in = 1.0 is the internally consistent
#: choice; it also matches the AmberTools MMPBSA.py default. A value of 2.0 was
#: used before 2026-08-18, which damped the polar term by roughly half relative
#: to the Coulomb term it is meant to balance.
GB_SOLUTE_DIELECTRIC = 1.0


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
    gb.setSoluteDielectric(float(GB_SOLUTE_DIELECTRIC))
    if include_sa:
        gb.setSurfaceAreaEnergy(2.25936)  # kJ/mol/nm^2
    else:
        gb.setSurfaceAreaEnergy(0.0)

    scales = _obc_scale_factors()
    fallback = scales["_default"]
    for atom in topology.atoms():
        q, _sigma, _eps = nb.getParticleParameters(atom.index)
        symbol = atom.element.symbol if atom.element is not None else atom.name[0]
        gb.addParticle(
            q,
            _radius_from_symbol(symbol),
            scales.get(symbol.upper(), fallback),
        )

    gb.setForceGroup(2)
    system.addForce(gb)
    return system


def _make_context(system):
    openmm = require_module("openmm")
    unit = require_module("openmm.unit")
    integrator = openmm.VerletIntegrator(0.001 * unit.picoseconds)
    platform_name = os.environ.get("OPENMM_PLATFORM", "CPU").strip() or "CPU"
    try:
        platform = openmm.Platform.getPlatformByName(platform_name)
        properties = {}
        if platform_name == "CPU":
            properties["Threads"] = os.environ.get("OPENMM_CPU_THREADS", "1")
        return openmm.Context(system, integrator, platform, properties)
    except Exception:
        return openmm.Context(system, integrator)


def _build_h_relax_context(topology, forcefield, k_kj_per_nm2: float = 10_000.0):
    """Build a reusable OpenMM context for H-atom relaxation.

    Returns (context, restraint_force, heavy_particle_list) where
    heavy_particle_list is a list of (atom_index, restraint_particle_index)
    pairs.  Call _apply_h_relax to minimise a specific snapshot without
    rebuilding the context.
    """
    app = require_module("openmm.app")
    openmm = require_module("openmm")
    unit = require_module("openmm.unit")

    system = forcefield.createSystem(
        topology,
        nonbondedMethod=app.NoCutoff,
        constraints=None,
    )

    restraint = openmm.CustomExternalForce("k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
    restraint.addGlobalParameter("k", k_kj_per_nm2)
    restraint.addPerParticleParameter("x0")
    restraint.addPerParticleParameter("y0")
    restraint.addPerParticleParameter("z0")
    heavy_particles: list[tuple[int, int]] = []
    for atom in topology.atoms():
        if atom.element is None or atom.element.symbol.upper() != "H":
            rp_idx = restraint.addParticle(atom.index, [0.0, 0.0, 0.0])
            heavy_particles.append((atom.index, rp_idx))
    system.addForce(restraint)

    integrator = openmm.VerletIntegrator(0.001 * unit.picoseconds)
    platform_name = os.environ.get("OPENMM_PLATFORM", "CPU").strip() or "CPU"
    try:
        platform = openmm.Platform.getPlatformByName(platform_name)
        properties = {}
        if platform_name == "CPU":
            properties["Threads"] = os.environ.get("OPENMM_CPU_THREADS", "1")
        context = openmm.Context(system, integrator, platform, properties)
    except Exception:
        context = openmm.Context(system, integrator)
    return context, restraint, heavy_particles


def _apply_h_relax(
    context,
    restraint,
    heavy_particles: list[tuple[int, int]],
    positions_nm: np.ndarray,
    max_iters: int = 200,
    tolerance_kj_per_nm: float = 100.0,
) -> np.ndarray:
    """Apply H-atom relaxation using a pre-built context (fast, reusable).

    Updates the restraint reference positions for this snapshot, minimises,
    and returns the corrected coordinates.  Only H atoms move; all non-H
    atoms are harmonically restrained to their MD positions.
    """
    openmm = require_module("openmm")
    unit = require_module("openmm.unit")

    # Update per-particle restraint reference positions for this snapshot.
    for atom_idx, rp_idx in heavy_particles:
        x, y, z = positions_nm[atom_idx]
        restraint.setParticleParameters(rp_idx, atom_idx, [float(x), float(y), float(z)])
    restraint.updateParametersInContext(context)

    pos = [
        openmm.Vec3(float(x), float(y), float(z)) * unit.nanometer
        for x, y, z in positions_nm
    ]
    context.setPositions(pos)
    openmm.LocalEnergyMinimizer.minimize(
        context, tolerance=tolerance_kj_per_nm, maxIterations=max_iters
    )
    state = context.getState(getPositions=True)
    return np.array(state.getPositions().value_in_unit(unit.nanometer))


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
    sample_window_ns: float | None = None,
    total_time_ns: float | None = None,
    sample_last_frames: int | None = None,
    allowed_frames: Iterable[int] | None = None,
) -> np.ndarray:
    if allowed_frames is not None:
        allowed = np.asarray(sorted({int(f) for f in allowed_frames if 0 <= int(f) < n_frames}), dtype=int)
        if allowed.size == 0:
            raise ValueError(
                "No usable frames remain after screening; widen the sampling window "
                "or relax the contact threshold."
            )
        # Take the most recent frames, matching the terminal-window convention
        # used by the unscreened paths.
        take = min(max(1, int(n_snapshots)), allowed.size)
        if sample_last_frames is not None and int(sample_last_frames) > 0:
            take = min(take, int(sample_last_frames))
        return allowed[-take:]

    if sample_last_frames is not None and int(sample_last_frames) > 0:
        start = max(0, int(n_frames) - int(sample_last_frames))
        available = n_frames - start
        n_take = min(max(1, int(n_snapshots)), max(1, available))
        return np.linspace(start, n_frames - 1, num=n_take, dtype=int)

    if (
        sample_window_ns is not None
        and sample_window_ns > 0.0
        and total_time_ns is not None
        and np.isfinite(total_time_ns)
        and float(total_time_ns) > 0.0
        and n_frames > 1
    ):
        window_ns = float(sample_window_ns)
        total_ns = float(total_time_ns)
        if window_ns >= total_ns:
            start = 0
        else:
            frac_start = max(0.0, (total_ns - window_ns) / total_ns)
            start = int(np.ceil(frac_start * float(n_frames - 1)))
        start = min(max(0, start), max(0, n_frames - 1))
        available = n_frames - start
        n_take = min(max(1, int(n_snapshots)), max(1, available))
        return np.linspace(start, n_frames - 1, num=n_take, dtype=int)

    # Some legacy DCDs carry inflated dt metadata due to an interval double-count.
    # Example symptom: ~50000 ps/frame when the actual spacing is ~50 ps/frame.
    max_reasonable_dt_ps = 1_000.0
    corrected_dt_ps = dt_ps
    if corrected_dt_ps is not None and np.isfinite(corrected_dt_ps) and corrected_dt_ps > max_reasonable_dt_ps:
        candidate = float(corrected_dt_ps) / 1000.0
        if candidate <= max_reasonable_dt_ps:
            corrected_dt_ps = candidate
        else:
            corrected_dt_ps = None

    if (
        sample_window_ns is not None
        and sample_window_ns > 0.0
        and corrected_dt_ps is not None
        and np.isfinite(corrected_dt_ps)
        and corrected_dt_ps > 0.0
        and corrected_dt_ps <= max_reasonable_dt_ps
        and n_frames > 1
    ):
        total_time_ps = float(n_frames - 1) * float(corrected_dt_ps)
        window_ps = float(sample_window_ns) * 1000.0
        start_time_ps = max(0.0, total_time_ps - window_ps)
        start = int(np.ceil(start_time_ps / float(corrected_dt_ps)))
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
    sample_window_ns: float | None = None,
    total_time_ns: float | None = None,
    sample_last_frames: int | None = None,
    analysis_topology_pdb_path: Path | None = None,
    allowed_frames: Iterable[int] | None = None,
) -> MMGBSAResult:
    """Compute MM/GBSA-style decomposition from explicit-MD snapshots.

    Note: If analysis_topology_pdb_path is provided (solute-only topology),
    it will be used for BOTH force field setup AND trajectory loading.
    The minimized_pdb_path is ignored in this case.

    ``allowed_frames`` restricts sampling to a whitelist of frame indices --
    used to exclude frames carrying unphysical ligand-protein contacts, which
    would otherwise contribute a huge spurious repulsive term. When given it
    overrides the window-based selection and the most recent allowed frames are
    taken.
    """
    import MDAnalysis as mda

    app = require_module("openmm.app")

    # Use analysis topology if provided (solute-only), otherwise minimized PDB
    topology_pdb_path = analysis_topology_pdb_path if analysis_topology_pdb_path is not None else minimized_pdb_path

    with open(topology_pdb_path, "r") as handle:
        solute_pdb = app.PDBFile(handle)

    ligand = load_ligand_molecule(ligand_sdf)
    forcefield = build_forcefield([ligand])

    complex_top = solute_pdb.topology

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

    # Build a single reusable H-relax context (avoids rebuilding per snapshot).
    h_relax_ctx, h_relax_force, h_relax_heavy = _build_h_relax_context(complex_top, forcefield)

    # Load trajectory with the same topology used for force field setup
    u = mda.Universe(str(topology_pdb_path), str(trajectory_dcd_path))
    n_frames = len(u.trajectory)
    if n_frames < 1:
        raise ValueError("Empty trajectory for MM/GBSA evaluation.")

    complex_n = u.atoms.n_atoms

    snap_idx = _select_snapshot_indices(
        n_frames=n_frames,
        discard_fraction=discard_fraction,
        n_snapshots=n_snapshots,
        dt_ps=getattr(u.trajectory, "dt", None),
        sample_window_ns=sample_window_ns,
        total_time_ns=total_time_ns,
        sample_last_frames=sample_last_frames,
        allowed_frames=allowed_frames,
    )

    d_vdw: list[float] = []
    d_elec: list[float] = []
    d_gb: list[float] = []
    d_sa: list[float] = []
    d_tot: list[float] = []
    absolute: dict[str, list[float]] = {
        f"{phase}_{term}": []
        for phase in ("complex", "receptor", "ligand")
        for term in ("vdw", "elec", "gb", "sa", "total")
    }

    for frame in snap_idx:
        u.trajectory[int(frame)]
        pos_a = np.asarray(u.atoms.positions, dtype=float)
        if pos_a.shape[0] != complex_n:
            raise ValueError(
                f"Trajectory has {pos_a.shape[0]} atoms but expected {complex_n}."
            )
        solute_nm = pos_a / 10.0
        # Relax H-atom positions to remove finite-timestep Langevin artifacts
        # (e.g., sub-Å inter-molecular H–H overlaps from SHAKE-constrained MD).
        # Only H atoms move; heavy atoms are restrained to their MD coordinates.
        solute_nm = _apply_h_relax(h_relax_ctx, h_relax_force, h_relax_heavy, solute_nm)
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

        for phase, (ev, ee, gp, gt) in (
            ("complex", (e_vdw_c, e_elec_c, gb_p_c, gb_t_c)),
            ("receptor", (e_vdw_r, e_elec_r, gb_p_r, gb_t_r)),
            ("ligand", (e_vdw_l, e_elec_l, gb_p_l, gb_t_l)),
        ):
            sa_term = gt - gp
            absolute[f"{phase}_vdw"].append(float(ev))
            absolute[f"{phase}_elec"].append(float(ee))
            absolute[f"{phase}_gb"].append(float(gp))
            absolute[f"{phase}_sa"].append(float(sa_term))
            absolute[f"{phase}_total"].append(float(ev + ee + gp + sa_term))

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

    # Explicitly clean up OpenMM contexts to prevent memory leaks and segfaults
    del h_relax_ctx
    for ctx in contexts.values():
        del ctx
    contexts.clear()
    for sys in systems.values():
        del sys
    systems.clear()

    absolute_summary = tuple((key, _stats(values)) for key, values in absolute.items())

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
        absolute_terms=absolute_summary,
    )
