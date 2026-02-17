from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor, as_completed
import json
import logging
import os
from pathlib import Path
import re
import time
import warnings

import numpy as np
import pandas as pd
from scipy import stats

from ..md.manifest import load_manifest

_PROJECT_ROOT = Path(__file__).resolve().parents[2]


def _format_seconds(seconds: float) -> str:
    seconds = max(0.0, float(seconds))
    mins, sec = divmod(int(seconds), 60)
    hrs, mins = divmod(mins, 60)
    if hrs > 0:
        return f"{hrs:d}h {mins:02d}m {sec:02d}s"
    return f"{mins:02d}m {sec:02d}s"


def _default_profile_workers() -> int:
    cpu = os.cpu_count() or 2
    return max(1, min(8, cpu - 1))


def _silence_mdanalysis_noise() -> None:
    logging.getLogger("MDAnalysis").setLevel(logging.WARNING)
    logging.getLogger("MDAnalysis.topology.guessers").setLevel(logging.WARNING)
    logging.getLogger("MDAnalysis.topology.PDBParser").setLevel(logging.WARNING)
    warnings.filterwarnings(
        "ignore",
        message="DCDReader currently makes independent timesteps*",
        category=DeprecationWarning,
    )



def _nonempty_path(value: object) -> Path | None:
    text = str(value or "").strip()
    if not text:
        return None
    return Path(text)


def _remap_to_local_workspace(candidate: Path | None) -> Path | None:
    """Map stale absolute paths (e.g., /scratch/.../nnrti-mechanisms/...) to local repo."""
    if candidate is None:
        return None
    if candidate.exists():
        return candidate

    marker = "nnrti-mechanisms/"
    text = str(candidate)
    if marker not in text:
        return candidate

    rel = text.split(marker, 1)[1]
    mapped = _PROJECT_ROOT / rel
    if mapped.exists():
        return mapped
    return candidate


def _resolve_local_path(candidate: Path | None, fallback: Path | None = None) -> Path | None:
    candidate = _remap_to_local_workspace(candidate)
    fallback = _remap_to_local_workspace(fallback)
    if candidate is not None and candidate.exists():
        return candidate
    if fallback is not None and fallback.exists():
        return fallback
    return candidate or fallback


def _infer_rep_dir(row: pd.Series) -> Path:
    for key in ("analysis_dcd", "trajectory_dcd", "prepared_topology_pdb", "minimized_pdb"):
        val = str(row.get(key) or "").strip()
        if val:
            return Path(val).parent
    return Path(".")


def _prepare_profile_jobs(run_df: pd.DataFrame) -> list[dict]:
    jobs: list[dict] = []
    for _, row in run_df.iterrows():
        rep_dir = _infer_rep_dir(row)
        safe = str(row["safe_label"])
        rep = int(row["replicate"])
        topo = _resolve_local_path(
            _nonempty_path(row.get("analysis_topology_pdb")),
            rep_dir / f"{safe}_rep{rep:02d}_analysis_topology.pdb",
        )
        dcd = _resolve_local_path(
            _nonempty_path(row.get("analysis_dcd")),
            rep_dir / f"{safe}_rep{rep:02d}_analysis.dcd",
        )
        if topo is None or dcd is None or not topo.exists() or not dcd.exists():
            continue
        jobs.append(
            {
                "structure": str(row["structure"]),
                "mutation": str(row["mutation"]),
                "safe_label": safe,
                "replicate": rep,
                "topology": str(topo),
                "trajectory": str(dcd),
            }
        )
    return jobs


def _run_profile_jobs(jobs: list[dict], worker_fn, label: str, workers: int | None = None) -> list[dict]:
    if not jobs:
        return []
    n_workers = _default_profile_workers() if workers is None else max(1, int(workers))
    total = len(jobs)
    start = time.time()
    rows: list[dict] = []
    logging.info("%s: %d trajectory jobs, workers=%d", label, total, n_workers)

    def _log_progress(done: int) -> None:
        if done == 1 or done == total or done % max(1, total // 10) == 0:
            elapsed = time.time() - start
            rate = done / elapsed if elapsed > 0 else 0.0
            eta = ((total - done) / rate) if rate > 0 else float("inf")
            eta_txt = _format_seconds(eta) if np.isfinite(eta) else "unknown"
            logging.info(
                "%s progress %d/%d (%.1f%%) elapsed=%s eta=%s",
                label,
                done,
                total,
                100.0 * done / total,
                _format_seconds(elapsed),
                eta_txt,
            )

    if n_workers == 1:
        done = 0
        for job in jobs:
            out_rows, err = worker_fn(job)
            rows.extend(out_rows)
            if err is not None:
                logging.warning(
                    "%s failed for %s rep%d: %s",
                    label,
                    job["mutation"],
                    int(job["replicate"]),
                    err,
                )
            done += 1
            _log_progress(done)
        return rows

    with ThreadPoolExecutor(max_workers=n_workers) as pool:
        futures = {pool.submit(worker_fn, job): job for job in jobs}
        done = 0
        for fut in as_completed(futures):
            job = futures[fut]
            try:
                out_rows, err = fut.result()
            except Exception as exc:
                out_rows, err = [], str(exc)
            rows.extend(out_rows)
            if err is not None:
                logging.warning(
                    "%s failed for %s rep%d: %s",
                    label,
                    job["mutation"],
                    int(job["replicate"]),
                    err,
                )
            done += 1
            _log_progress(done)
    return rows


def _ca_rmsd_worker(job: dict, frame_stride: int, max_frames: int) -> tuple[list[dict], str | None]:
    _silence_mdanalysis_noise()
    try:
        import MDAnalysis as mda
        from MDAnalysis.analysis import align
    except Exception as exc:
        return [], str(exc)

    try:
        u = mda.Universe(job["topology"], job["trajectory"])
        ref = mda.Universe(job["topology"])
        prot_all = u.select_atoms("protein")
        try:
            from MDAnalysis import transformations as trans

            u.trajectory.add_transformations(
                trans.NoJump(check_continuity=False),
                trans.center_in_box(prot_all, center="geometry", wrap=False),
            )
        except Exception:
            pass

        out: list[dict] = []
        kept = 0
        stride = max(1, int(frame_stride))
        for idx, _ in enumerate(u.trajectory):
            if idx % stride != 0:
                continue
            if kept >= int(max_frames):
                break
            align.alignto(u, ref, select="protein and name CA", weights="mass")
            ca = u.select_atoms("protein and name CA")
            ca_ref = ref.select_atoms("protein and name CA")
            if ca.n_atoms == 0 or ca_ref.n_atoms == 0:
                break
            diff = ca.positions - ca_ref.positions
            rmsd = float(np.sqrt(np.mean(np.sum(diff * diff, axis=1))))
            out.append(
                {
                    "structure": job["structure"],
                    "mutation": job["mutation"],
                    "safe_label": job["safe_label"],
                    "replicate": int(job["replicate"]),
                    "frame_index": int(u.trajectory.frame),
                    "time_ps": float(getattr(u.trajectory.ts, "time", np.nan)),
                    "ca_rmsd_angstrom": rmsd,
                }
            )
            kept += 1
        return out, None
    except Exception as exc:
        return [], str(exc)


def _com_distance_worker(
    job: dict,
    ligand_resname: str,
    frame_stride: int,
    max_frames: int,
) -> tuple[list[dict], str | None]:
    _silence_mdanalysis_noise()
    try:
        import MDAnalysis as mda
        from MDAnalysis.lib.distances import distance_array
    except Exception as exc:
        return [], str(exc)

    try:
        u = mda.Universe(job["topology"], job["trajectory"])
        lig = u.select_atoms(f"resname {ligand_resname}")
        prot = u.select_atoms(f"protein and not resname {ligand_resname}")
        if lig.n_atoms == 0 or prot.n_atoms == 0:
            return [], f"empty selection (lig={lig.n_atoms}, prot={prot.n_atoms})"
        try:
            from MDAnalysis import transformations as trans

            u.trajectory.add_transformations(
                trans.NoJump(check_continuity=False),
                trans.center_in_box(prot, center="geometry", wrap=False),
            )
        except Exception:
            pass

        out: list[dict] = []
        kept = 0
        stride = max(1, int(frame_stride))
        for idx, _ in enumerate(u.trajectory):
            if idx % stride != 0:
                continue
            if kept >= int(max_frames):
                break
            lig_com = np.asarray(lig.center_of_mass(), dtype=float).reshape(1, 3)
            prot_com = np.asarray(prot.center_of_mass(), dtype=float).reshape(1, 3)
            d = float(distance_array(lig_com, prot_com, box=u.dimensions).min())
            out.append(
                {
                    "structure": job["structure"],
                    "mutation": job["mutation"],
                    "safe_label": job["safe_label"],
                    "replicate": int(job["replicate"]),
                    "frame_index": int(u.trajectory.frame),
                    "time_ps": float(getattr(u.trajectory.ts, "time", np.nan)),
                    "com_distance_angstrom": d,
                }
            )
            kept += 1
        return out, None
    except Exception as exc:
        return [], str(exc)


def _pocket_volume_worker(
    job: dict,
    ligand_resname: str,
    frame_stride: int,
    max_frames: int,
    grid_spacing: float,
    pocket_radius_angstrom: float,
) -> tuple[list[dict], str | None]:
    _silence_mdanalysis_noise()
    try:
        import MDAnalysis as mda
        from .metrics import pocket_volume_proxy_from_universe
    except Exception as exc:
        return [], str(exc)

    try:
        u = mda.Universe(job["topology"], job["trajectory"])
        lig = u.select_atoms(f"resname {ligand_resname}")
        prot = u.select_atoms("protein")
        if lig.n_atoms == 0 or prot.n_atoms == 0:
            return [], f"empty selection (lig={lig.n_atoms}, prot={prot.n_atoms})"
        try:
            from MDAnalysis import transformations as trans

            u.trajectory.add_transformations(
                trans.NoJump(check_continuity=False),
                trans.center_in_box(prot, center="geometry", wrap=False),
            )
        except Exception:
            pass

        out: list[dict] = []
        kept = 0
        stride = max(1, int(frame_stride))
        for idx, _ in enumerate(u.trajectory):
            if idx % stride != 0:
                continue
            if kept >= int(max_frames):
                break
            v = pocket_volume_proxy_from_universe(
                u,
                ligand_resname=ligand_resname,
                grid_spacing=float(grid_spacing),
                radius_angstrom=float(pocket_radius_angstrom),
            )
            out.append(
                {
                    "structure": job["structure"],
                    "mutation": job["mutation"],
                    "safe_label": job["safe_label"],
                    "replicate": int(job["replicate"]),
                    "frame_index": int(u.trajectory.frame),
                    "time_ps": float(getattr(u.trajectory.ts, "time", np.nan)),
                    "pocket_volume_proxy_angstrom3": float(v),
                    "grid_spacing_angstrom": float(grid_spacing),
                    "pocket_radius_angstrom": float(pocket_radius_angstrom),
                }
            )
            kept += 1
        return out, None
    except Exception as exc:
        return [], str(exc)


def collect_md_results(manifest_path: Path, md_results_dir: Path | None = None) -> pd.DataFrame:
    """Collect per-replicate MD run metadata from JSON files."""
    del md_results_dir
    tasks = load_manifest(manifest_path)
    rows: list[dict] = []

    for task in tasks:
        json_path = Path(task.output_json)
        if not json_path.exists():
            logging.warning("Missing result for task %d: %s", task.task_id, json_path)
            continue

        try:
            data = json.loads(json_path.read_text())
        except Exception as exc:
            logging.error("Invalid JSON for task %d: %s", task.task_id, exc)
            continue

        rows.append(
            {
                "task_id": task.task_id,
                "structure": task.structure,
                "mutation": task.mutation,
                "safe_label": task.safe_label,
                "replicate": int(task.replicate),
                "fold_reduction": task.fold_reduction,
                "minimized_pdb": data.get("minimized_pdb") or task.minimized_pdb,
                "prepared_topology_pdb": data.get("prepared_topology_pdb") or task.prepared_topology_pdb,
                "prepared_system_xml": data.get("prepared_system_xml") or task.prepared_system_xml,
                "analysis_dcd": data.get("analysis_dcd", ""),
                "analysis_topology_pdb": data.get("analysis_topology_pdb", ""),
                "final_pdb": data.get("final_pdb", ""),
                "status": data.get("status", ""),
                "elapsed_seconds": data.get("elapsed_seconds"),
                "ligand_sdf": data.get("ligand_sdf") or task.ligand_sdf,
                "ligand_resname": data.get("ligand_resname") or task.ligand_resname,
            }
        )

    return pd.DataFrame(rows)


def collect_ca_rmsd_profiles(
    run_df: pd.DataFrame,
    frame_stride: int = 10,
    max_frames: int = 400,
    workers: int | None = None,
) -> pd.DataFrame:
    if run_df.empty:
        return pd.DataFrame()
    jobs = _prepare_profile_jobs(run_df)
    if not jobs:
        return pd.DataFrame()

    def _runner(job: dict) -> tuple[list[dict], str | None]:
        return _ca_rmsd_worker(job, frame_stride=frame_stride, max_frames=max_frames)

    rows = _run_profile_jobs(jobs, _runner, label="RMSD profiles", workers=workers)
    return pd.DataFrame(rows)


def collect_com_distance_profiles(
    run_df: pd.DataFrame,
    ligand_resname: str,
    frame_stride: int = 10,
    max_frames: int = 400,
    workers: int | None = None,
) -> pd.DataFrame:
    """Collect DOR-vs-protein center-of-mass distance profiles."""
    if run_df.empty:
        return pd.DataFrame()
    jobs = _prepare_profile_jobs(run_df)
    if not jobs:
        return pd.DataFrame()

    def _runner(job: dict) -> tuple[list[dict], str | None]:
        return _com_distance_worker(
            job,
            ligand_resname=ligand_resname,
            frame_stride=frame_stride,
            max_frames=max_frames,
        )

    rows = _run_profile_jobs(jobs, _runner, label="COM-distance profiles", workers=workers)
    return pd.DataFrame(rows)


def collect_pocket_volume_profiles(
    run_df: pd.DataFrame,
    ligand_resname: str,
    frame_stride: int = 5,
    max_frames: int = 400,
    grid_spacing: float = 0.75,
    pocket_radius_angstrom: float = 8.0,
    workers: int | None = None,
) -> pd.DataFrame:
    """Collect per-frame pocket-volume proxy traces."""
    if run_df.empty:
        return pd.DataFrame()
    jobs = _prepare_profile_jobs(run_df)
    if not jobs:
        return pd.DataFrame()

    def _runner(job: dict) -> tuple[list[dict], str | None]:
        return _pocket_volume_worker(
            job,
            ligand_resname=ligand_resname,
            frame_stride=frame_stride,
            max_frames=max_frames,
            grid_spacing=grid_spacing,
            pocket_radius_angstrom=pocket_radius_angstrom,
        )

    rows = _run_profile_jobs(jobs, _runner, label="Pocket-volume profiles", workers=workers)
    return pd.DataFrame(rows)


def compute_boundness_qc(
    run_df: pd.DataFrame,
    ligand_resname: str,
    traj_frame_stride: int = 20,
    traj_max_frames: int = 50,
    bound_threshold_angstrom: float = 6.0,
) -> pd.DataFrame:
    if run_df.empty:
        return pd.DataFrame()

    has_mda = False
    try:
        import MDAnalysis  # noqa: F401
        from MDAnalysis.lib.distances import capped_distance  # noqa: F401

        has_mda = True
    except Exception:
        has_mda = False

    rows: list[dict] = []
    for _, row in run_df.iterrows():
        rep_dir = _infer_rep_dir(row)
        safe = str(row["safe_label"])
        rep = int(row["replicate"])

        topo = _resolve_local_path(
            _nonempty_path(row.get("analysis_topology_pdb")),
            rep_dir / f"{safe}_rep{rep:02d}_analysis_topology.pdb",
        )
        dcd = _resolve_local_path(
            _nonempty_path(row.get("analysis_dcd")),
            rep_dir / f"{safe}_rep{rep:02d}_analysis.dcd",
        )

        start_min_dist = float("nan")
        traj_min_dist = float("nan")
        n_frames = 0

        if has_mda and topo is not None and topo.exists():
            try:
                import MDAnalysis
                from MDAnalysis.lib.distances import capped_distance

                # Compute start-structure min distance from topology PDB
                u_start = MDAnalysis.Universe(str(topo))
                lig_s = u_start.select_atoms(f"resname {ligand_resname} and not name H*")
                prot_s = u_start.select_atoms(f"protein and not resname {ligand_resname} and not name H*")
                if lig_s.n_atoms > 0 and prot_s.n_atoms > 0:
                    _pairs_s, dists_s = capped_distance(
                        lig_s.positions, prot_s.positions,
                        max_cutoff=30.0, box=u_start.dimensions, return_distances=True,
                    )
                    start_min_dist = float(np.min(dists_s)) if len(dists_s) else float("nan")

                # Compute trajectory min distance if DCD available
                if dcd is not None and dcd.exists():
                    u = MDAnalysis.Universe(str(topo), str(dcd))
                    lig = u.select_atoms(f"resname {ligand_resname} and not name H*")
                    prot = u.select_atoms(f"protein and not resname {ligand_resname} and not name H*")
                    mins: list[float] = []
                    for i, _ in enumerate(u.trajectory[:: max(1, traj_frame_stride)]):
                        if i >= traj_max_frames:
                            break
                        n_frames += 1
                        _pairs, dists = capped_distance(
                            lig.positions, prot.positions,
                            max_cutoff=30.0, box=u.dimensions, return_distances=True,
                        )
                        mins.append(float(np.min(dists)) if len(dists) else 30.0)
                    if mins:
                        traj_min_dist = float(np.min(mins))
            except Exception:
                pass

        bound_start = bool(np.isfinite(start_min_dist) and start_min_dist <= bound_threshold_angstrom)
        bound_traj = bool(np.isfinite(traj_min_dist) and traj_min_dist <= bound_threshold_angstrom)

        if not bound_start:
            qc_flag = "UNBOUND_START"
        elif np.isfinite(traj_min_dist) and not bound_traj:
            qc_flag = "UNBOUND_TRAJECTORY"
        else:
            qc_flag = "OK"

        rows.append(
            {
                "structure": row["structure"],
                "mutation": row["mutation"],
                "safe_label": safe,
                "replicate": rep,
                "min_distance_start_angstrom": start_min_dist,
                "min_distance_trajectory_angstrom": traj_min_dist,
                "trajectory_frames_sampled": n_frames,
                "bound_threshold_angstrom": float(bound_threshold_angstrom),
                "is_bound_start": bound_start,
                "is_bound_trajectory": bound_traj if np.isfinite(traj_min_dist) else np.nan,
                "qc_flag": qc_flag,
                "topology_pdb_used": str(topo) if topo is not None else "",
                "trajectory_dcd_used": str(dcd) if dcd is not None else "",
            }
        )

    return pd.DataFrame(rows)


def compute_structural_metrics(
    run_df: pd.DataFrame,
    ligand_resname: str,
    frame_stride: int = 5,
    max_frames: int = 200,
) -> pd.DataFrame:
    from .metrics import compute_ensemble_metrics

    if run_df.empty:
        return pd.DataFrame()

    def _infer_total_ns_from_state_csv(rep_dir: Path, safe_label: str, replicate: int) -> float | None:
        # Use the OpenMM StateDataReporter output to infer the true production length.
        state_csv = rep_dir / f"{safe_label}_rep{int(replicate):02d}_md_state.csv"
        if not state_csv.exists():
            return None
        try:
            sdf = pd.read_csv(state_csv)
        except Exception:
            return None
        step_col = None
        for c in ('#"Step"', "Step"):
            if c in sdf.columns:
                step_col = c
                break
        if step_col is None or sdf.empty:
            return None
        max_step = pd.to_numeric(sdf[step_col], errors="coerce").dropna()
        if max_step.empty:
            return None
        # 2 fs timestep => ns = steps * 2 fs / 1e6 fs/ns
        return float(max_step.max()) * 2.0 / 1_000_000.0

    rows: list[dict] = []
    for _, row in run_df.iterrows():
        rep_dir = _infer_rep_dir(row)
        safe = str(row["safe_label"])
        rep = int(row["replicate"])

        topo = _resolve_local_path(
            _nonempty_path(row.get("analysis_topology_pdb")),
            rep_dir / f"{safe}_rep{rep:02d}_analysis_topology.pdb",
        )
        dcd = _resolve_local_path(
            _nonempty_path(row.get("analysis_dcd")),
            rep_dir / f"{safe}_rep{rep:02d}_analysis.dcd",
        )
        if topo is None or dcd is None or not topo.exists() or not dcd.exists():
            logging.warning("Missing trajectory inputs for %s rep%d", row["mutation"], rep)
            continue

        total_ns = _infer_total_ns_from_state_csv(rep_dir, safe, rep)
        time_source = "md_state_csv" if total_ns is not None and np.isfinite(total_ns) else "trajectory_dt"

        try:
            ens = compute_ensemble_metrics(
                topology_pdb_path=topo,
                trajectory_dcd_path=dcd,
                ligand_resname=ligand_resname,
                frame_stride=frame_stride,
                max_frames=max_frames,
                total_time_ns=total_ns,
            )
            rows.append(
                {
                    "structure": row["structure"],
                    "mutation": row["mutation"],
                    "safe_label": safe,
                    "replicate": rep,
                    "contact_count": ens.contact_count_mean,
                    "contact_count_std": ens.contact_count_std,
                    "hbond_count": ens.hbond_count_mean,
                    "hbond_count_std": ens.hbond_count_std,
                    "pocket_volume_proxy": ens.pocket_volume_proxy_mean,
                    "pocket_volume_proxy_std": ens.pocket_volume_proxy_std,
                    "metric_n_frames": ens.n_frames,
                    "metric_source": "trajectory",
                    "metric_sample_window_ns": 1.0,
                    "metric_time_source": time_source,
                    "fold_reduction": row["fold_reduction"],
                }
            )
        except Exception as exc:
            logging.error("Structural metrics failed for %s rep%d: %s", row["mutation"], rep, exc)

    return pd.DataFrame(rows)


def compute_mmgbsa_metrics(
    run_df: pd.DataFrame,
    ligand_resname: str,
    n_snapshots: int = 100,
    discard_fraction: float = 0.25,
) -> pd.DataFrame:
    import gc
    from ..md.openmm.mmgbsa import compute_mmgbsa_from_trajectory

    if run_df.empty:
        return pd.DataFrame()

    rows: list[dict] = []
    for idx, row in run_df.iterrows():
        rep_dir = _infer_rep_dir(row)
        safe = str(row["safe_label"])
        rep = int(row["replicate"])

        min_pdb = _resolve_local_path(
            _nonempty_path(row.get("minimized_pdb")),
            rep_dir / f"{safe}_minimized_rep{rep:02d}.pdb",
        )
        dcd = _resolve_local_path(
            _nonempty_path(row.get("analysis_dcd")),
            rep_dir / f"{safe}_rep{rep:02d}_analysis.dcd",
        )
        analysis_topo = _resolve_local_path(
            _nonempty_path(row.get("analysis_topology_pdb")),
            rep_dir / f"{safe}_rep{rep:02d}_analysis_topology.pdb",
        )
        ligand_sdf = _resolve_local_path(_nonempty_path(row.get("ligand_sdf")))

        if None in {min_pdb, dcd, analysis_topo, ligand_sdf}:
            logging.warning("Missing MM/GBSA inputs for %s rep%d", row["mutation"], rep)
            continue
        if not min_pdb.exists() or not dcd.exists() or not analysis_topo.exists() or not ligand_sdf.exists():
            logging.warning("Unavailable MM/GBSA paths for %s rep%d", row["mutation"], rep)
            continue

        try:
            mm = compute_mmgbsa_from_trajectory(
                minimized_pdb_path=min_pdb,
                trajectory_dcd_path=dcd,
                ligand_resname=ligand_resname,
                ligand_sdf=ligand_sdf,
                n_snapshots=n_snapshots,
                discard_fraction=discard_fraction,
                analysis_topology_pdb_path=analysis_topo,
            )
            rows.append(
                {
                    "structure": row["structure"],
                    "mutation": row["mutation"],
                    "safe_label": safe,
                    "replicate": rep,
                    "fold_reduction": row["fold_reduction"],
                    "binding_dg": mm.binding_dg_mean,
                    "binding_dg_std": mm.binding_dg_std,
                    "binding_dg_sem": mm.binding_dg_sem,
                    "binding_dg_vdw": mm.delta_e_vdw_mean,
                    "binding_dg_vdw_std": mm.delta_e_vdw_std,
                    "binding_dg_vdw_sem": mm.delta_e_vdw_sem,
                    "binding_dg_electrostatic": mm.delta_e_elec_mean,
                    "binding_dg_electrostatic_std": mm.delta_e_elec_std,
                    "binding_dg_electrostatic_sem": mm.delta_e_elec_sem,
                    "binding_dg_gb": mm.delta_g_gb_mean,
                    "binding_dg_gb_std": mm.delta_g_gb_std,
                    "binding_dg_gb_sem": mm.delta_g_gb_sem,
                    "binding_dg_sa": mm.delta_g_sa_mean,
                    "binding_dg_sa_std": mm.delta_g_sa_std,
                    "binding_dg_sa_sem": mm.delta_g_sa_sem,
                    "mmgbsa_snapshots": mm.n_snapshots,
                }
            )
        except Exception as exc:
            logging.error("MM/GBSA failed for %s rep%d: %s", row["mutation"], rep, exc)

        # Force garbage collection after each replicate to prevent memory buildup
        gc.collect()

    return pd.DataFrame(rows)


def compute_binding_ddg(mmgbsa_df: pd.DataFrame) -> pd.DataFrame:
    if mmgbsa_df.empty:
        return pd.DataFrame()

    df = mmgbsa_df.copy()
    wt = df[df["mutation"] == "WT"].set_index(["structure", "replicate"])
    if wt.empty:
        df["wt_binding_dg"] = np.nan
        df["ddg"] = np.nan
        return df

    for col in ["binding_dg", "binding_dg_vdw", "binding_dg_electrostatic", "binding_dg_gb", "binding_dg_sa"]:
        wt_col = f"wt_{col}"
        ddg_col = col.replace("binding_dg", "ddg")
        lookup = wt[col] if col in wt.columns else pd.Series(dtype=float)
        df[wt_col] = df.apply(lambda r: lookup.get((r["structure"], r["replicate"]), np.nan), axis=1)
        df[ddg_col] = df[col] - df[wt_col]

    return df


def merge_with_structural_metrics(ddg_df: pd.DataFrame, structural_metrics_df: pd.DataFrame) -> pd.DataFrame:
    if structural_metrics_df.empty:
        return ddg_df
    merged = ddg_df.merge(
        structural_metrics_df,
        on=["structure", "mutation", "safe_label", "replicate"],
        how="left",
        suffixes=("", "_struct"),
    )
    if "fold_reduction_struct" in merged.columns:
        merged = merged.drop(columns=["fold_reduction_struct"])
    return merged


def _corr_rows(metric_name: str, x: pd.Series, y: pd.Series) -> list[dict]:
    mask = np.isfinite(x.values) & np.isfinite(y.values)
    if mask.sum() < 3:
        return []
    x_valid = x[mask]
    y_valid = y[mask]
    pearson_r, pearson_p = stats.pearsonr(x_valid, y_valid)
    spearman_rho, spearman_p = stats.spearmanr(x_valid, y_valid)
    return [
        {
            "metric": metric_name,
            "pearson_r": float(pearson_r),
            "pearson_r2": float(pearson_r**2),
            "pearson_pvalue": float(pearson_p),
            "spearman_rho": float(spearman_rho),
            "spearman_pvalue": float(spearman_p),
            "n_mutations": int(mask.sum()),
        }
    ]


def compute_correlations(ddg_df: pd.DataFrame) -> pd.DataFrame:
    mut_df = ddg_df[ddg_df["mutation"] != "WT"].copy()
    if mut_df.empty:
        return pd.DataFrame()

    wt_df = ddg_df[ddg_df["mutation"] == "WT"].set_index(["structure", "replicate"])
    for metric in ["contact_count", "hbond_count", "pocket_volume_proxy", "binding_dg"]:
        if metric in mut_df.columns and metric in wt_df.columns:
            lookup = wt_df[metric]
            mut_df[f"{metric}_delta"] = mut_df.apply(
                lambda r: r[metric] - lookup.get((r["structure"], r["replicate"]), np.nan),
                axis=1,
            )

    agg_cols: dict[str, tuple[str, str]] = {
        "fold_reduction": ("fold_reduction", "first"),
        "n_replicates": ("replicate", "nunique"),
    }
    for metric in ["ddg", "contact_count_delta", "hbond_count_delta", "pocket_volume_proxy_delta", "binding_dg_delta"]:
        if metric in mut_df.columns:
            agg_cols[f"{metric}_mean"] = (metric, "mean")

    by_mut = mut_df.groupby("mutation", as_index=False).agg(**agg_cols)
    if by_mut.empty:
        return pd.DataFrame()

    y = by_mut["fold_reduction"]
    rows: list[dict] = []
    for col in by_mut.columns:
        if col.endswith("_mean") and col not in {"fold_reduction_mean"}:
            rows.extend(_corr_rows(col.replace("_mean", ""), by_mut[col], y))
    return pd.DataFrame(rows)


def run_result_collection(
    manifest_path: Path,
    md_results_dir: Path | None = None,
    output_dir: Path | None = None,
    ligand_resname: str = "2KW",
    compute_structural: bool = True,
    metric_frame_stride: int = 5,
    metric_max_frames: int = 200,
    mmgbsa_snapshots: int = 100,
    mmgbsa_discard_fraction: float = 0.25,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    if output_dir is None:
        output_dir = manifest_path.parent
    output_dir.mkdir(parents=True, exist_ok=True)

    logging.info("Collecting MD run metadata from %s", md_results_dir or manifest_path.parent)
    run_df = collect_md_results(manifest_path, md_results_dir)
    if run_df.empty:
        raise ValueError("No completed MD results found")

    force_recompute = str(os.environ.get("NNRTI_FORCE_RECOMPUTE", "0")).strip().lower() in {
        "1",
        "true",
        "yes",
        "y",
    }
    expected_keys = run_df[["structure", "mutation", "safe_label", "replicate"]].drop_duplicates().copy()

    mmgbsa_cache = output_dir / "mmgbsa_replicate_metrics.csv"
    mmgbsa_df = pd.DataFrame()
    if not force_recompute and mmgbsa_cache.exists():
        try:
            cached = pd.read_csv(mmgbsa_cache)
            required = {"structure", "mutation", "safe_label", "replicate", "binding_dg", "fold_reduction"}
            if required.issubset(set(cached.columns)):
                cached["replicate"] = pd.to_numeric(cached["replicate"], errors="coerce").astype("Int64")
                cached = cached.dropna(subset=["replicate"]).copy()
                cached["replicate"] = cached["replicate"].astype(int)
                cached = cached.merge(
                    expected_keys,
                    on=["structure", "mutation", "safe_label", "replicate"],
                    how="inner",
                )
                if len(cached) == len(expected_keys):
                    mmgbsa_df = cached
                    logging.info("Reusing cached MM/GBSA metrics from %s", mmgbsa_cache)
                else:
                    logging.info(
                        "Cached MM/GBSA metrics incomplete (%d/%d rows); recomputing.",
                        len(cached),
                        len(expected_keys),
                    )
        except Exception as exc:
            logging.info("Could not read cached MM/GBSA metrics (%s); recomputing.", exc)

    if mmgbsa_df.empty:
        logging.info("Running MM/GBSA snapshot analysis")
        mmgbsa_df = compute_mmgbsa_metrics(
            run_df,
            ligand_resname=ligand_resname,
            n_snapshots=mmgbsa_snapshots,
            discard_fraction=mmgbsa_discard_fraction,
        )
        if mmgbsa_df.empty:
            raise ValueError("No MM/GBSA metrics could be computed")
        mmgbsa_df.to_csv(mmgbsa_cache, index=False)

    ddg_df = compute_binding_ddg(mmgbsa_df)

    struct_df = pd.DataFrame()
    if compute_structural:
        struct_cache = output_dir / "structural_metrics.csv"
        if not force_recompute and struct_cache.exists():
            try:
                cached = pd.read_csv(struct_cache)
                required = {
                    "structure",
                    "mutation",
                    "safe_label",
                    "replicate",
                    "contact_count",
                    "hbond_count",
                    "pocket_volume_proxy",
                    "metric_sample_window_ns",
                }
                if required.issubset(set(cached.columns)):
                    cached["replicate"] = pd.to_numeric(cached["replicate"], errors="coerce").astype("Int64")
                    cached = cached.dropna(subset=["replicate"]).copy()
                    cached["replicate"] = cached["replicate"].astype(int)
                    cached = cached.merge(
                        expected_keys,
                        on=["structure", "mutation", "safe_label", "replicate"],
                        how="inner",
                    )
                    # Defensive: drop duplicated keys if any (see NOTE below).
                    cached = (
                        cached.drop_duplicates(
                            subset=["structure", "mutation", "safe_label", "replicate"],
                            keep="first",
                        )
                        .reset_index(drop=True)
                    )
                    if len(cached) == len(expected_keys):
                        struct_df = cached
                        logging.info("Reusing cached structural metrics from %s", struct_cache)
                    else:
                        logging.info(
                            "Cached structural metrics incomplete (%d/%d rows); recomputing.",
                            len(cached),
                            len(expected_keys),
                        )
            except Exception as exc:
                logging.info("Could not read cached structural metrics (%s); recomputing.", exc)

        if struct_df.empty:
            logging.info("Computing ensemble structural metrics")
            struct_df = compute_structural_metrics(
                run_df,
                ligand_resname=ligand_resname,
                frame_stride=metric_frame_stride,
                max_frames=metric_max_frames,
            )
        if not struct_df.empty:
            # Defensive: a duplicated (mutation, replicate) row can silently propagate
            # into ddg_full.csv via merge multiplication. Keep the first occurrence.
            before = len(struct_df)
            struct_df = (
                struct_df.drop_duplicates(
                    subset=["structure", "mutation", "safe_label", "replicate"],
                    keep="first",
                )
                .reset_index(drop=True)
            )
            dropped = before - len(struct_df)
            if dropped:
                logging.warning(
                    "Dropped %d duplicate structural-metric rows based on (structure, mutation, safe_label, replicate).",
                    dropped,
                )
            struct_df.to_csv(struct_cache, index=False)
            ddg_df = merge_with_structural_metrics(ddg_df, struct_df)

    logging.info("Computing boundness QC")
    qc_cache = output_dir / "boundness_qc.csv"
    qc_df = pd.DataFrame()
    if not force_recompute and qc_cache.exists():
        try:
            qc_df = pd.read_csv(qc_cache)
            if not qc_df.empty:
                logging.info("Reusing cached boundness QC from %s", qc_cache)
        except Exception:
            qc_df = pd.DataFrame()
    if qc_df.empty:
        qc_df = compute_boundness_qc(run_df, ligand_resname)
        if not qc_df.empty:
            qc_df.to_csv(qc_cache, index=False)

    logging.info("Computing RMSD convergence profiles")
    rmsd_cache = output_dir / "rmsd_ca_profiles.csv"
    rmsd_df = pd.DataFrame()
    if not force_recompute and rmsd_cache.exists():
        try:
            rmsd_df = pd.read_csv(rmsd_cache)
            if not rmsd_df.empty:
                logging.info("Reusing cached RMSD profiles from %s", rmsd_cache)
        except Exception:
            rmsd_df = pd.DataFrame()
    if rmsd_df.empty:
        rmsd_df = collect_ca_rmsd_profiles(run_df)
        if not rmsd_df.empty:
            rmsd_df.to_csv(rmsd_cache, index=False)

    logging.info("Computing DOR-RT COM distance convergence profiles")
    com_cache = output_dir / "com_distance_profiles.csv"
    com_df = pd.DataFrame()
    if not force_recompute and com_cache.exists():
        try:
            com_df = pd.read_csv(com_cache)
            if not com_df.empty:
                logging.info("Reusing cached COM distance profiles from %s", com_cache)
        except Exception:
            com_df = pd.DataFrame()
    if com_df.empty:
        com_df = collect_com_distance_profiles(
            run_df,
            ligand_resname=ligand_resname,
            frame_stride=metric_frame_stride,
            max_frames=max(400, metric_max_frames),
        )
        if not com_df.empty:
            com_df.to_csv(com_cache, index=False)

    ddg_df.to_csv(output_dir / "ddg_full.csv", index=False)

    correlations = compute_correlations(ddg_df)
    correlations.to_csv(output_dir / "correlation_analysis.csv", index=False)

    return correlations, ddg_df
