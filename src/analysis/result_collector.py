from __future__ import annotations

import json
import logging
from pathlib import Path
import re

import numpy as np
import pandas as pd
from scipy import stats

from ..md.manifest import load_manifest

_PROJECT_ROOT = Path(__file__).resolve().parents[2]



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
) -> pd.DataFrame:
    try:
        import MDAnalysis as mda
        from MDAnalysis.analysis import align
    except Exception:
        return pd.DataFrame()

    if run_df.empty:
        return pd.DataFrame()

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
            continue

        try:
            u = mda.Universe(str(topo), str(dcd))
            ref = mda.Universe(str(topo))
            for i, _ in enumerate(u.trajectory[:: max(1, frame_stride)]):
                if i >= max_frames:
                    break
                align.alignto(u, ref, select="protein and name CA", weights="mass")
                ca = u.select_atoms("protein and name CA")
                ca_ref = ref.select_atoms("protein and name CA")
                if ca.n_atoms == 0 or ca_ref.n_atoms == 0:
                    break
                diff = ca.positions - ca_ref.positions
                rmsd = float(np.sqrt(np.mean(np.sum(diff * diff, axis=1))))
                rows.append(
                    {
                        "structure": row["structure"],
                        "mutation": row["mutation"],
                        "safe_label": safe,
                        "replicate": rep,
                        "frame_index": int(u.trajectory.frame),
                        "time_ps": float(getattr(u.trajectory.ts, "time", np.nan)),
                        "ca_rmsd_angstrom": rmsd,
                    }
                )
        except Exception as exc:
            logging.warning("RMSD profile failed for %s rep%d: %s", row["mutation"], rep, exc)

    return pd.DataFrame(rows)


def collect_com_distance_profiles(
    run_df: pd.DataFrame,
    ligand_resname: str,
    frame_stride: int = 10,
    max_frames: int = 400,
) -> pd.DataFrame:
    """Collect DOR-vs-protein center-of-mass distance profiles."""
    try:
        import MDAnalysis as mda
    except Exception:
        return pd.DataFrame()

    if run_df.empty:
        return pd.DataFrame()

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
            continue

        try:
            u = mda.Universe(str(topo), str(dcd))
            lig = u.select_atoms(f"resname {ligand_resname}")
            prot = u.select_atoms(f"protein and not resname {ligand_resname}")
            if lig.n_atoms == 0 or prot.n_atoms == 0:
                continue
            for i, _ in enumerate(u.trajectory[:: max(1, frame_stride)]):
                if i >= max_frames:
                    break
                d = float(np.linalg.norm(lig.center_of_mass() - prot.center_of_mass()))
                rows.append(
                    {
                        "structure": row["structure"],
                        "mutation": row["mutation"],
                        "safe_label": safe,
                        "replicate": rep,
                        "frame_index": int(u.trajectory.frame),
                        "time_ps": float(getattr(u.trajectory.ts, "time", np.nan)),
                        "com_distance_angstrom": d,
                    }
                )
        except Exception as exc:
            logging.warning("COM-distance profile failed for %s rep%d: %s", row["mutation"], rep, exc)

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

        try:
            ens = compute_ensemble_metrics(
                topology_pdb_path=topo,
                trajectory_dcd_path=dcd,
                ligand_resname=ligand_resname,
                frame_stride=frame_stride,
                max_frames=max_frames,
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

    logging.info("Running MM/GBSA snapshot analysis")
    mmgbsa_df = compute_mmgbsa_metrics(
        run_df,
        ligand_resname=ligand_resname,
        n_snapshots=mmgbsa_snapshots,
        discard_fraction=mmgbsa_discard_fraction,
    )
    if mmgbsa_df.empty:
        raise ValueError("No MM/GBSA metrics could be computed")
    mmgbsa_df.to_csv(output_dir / "mmgbsa_replicate_metrics.csv", index=False)

    ddg_df = compute_binding_ddg(mmgbsa_df)

    struct_df = pd.DataFrame()
    if compute_structural:
        logging.info("Computing ensemble structural metrics")
        struct_df = compute_structural_metrics(
            run_df,
            ligand_resname=ligand_resname,
            frame_stride=metric_frame_stride,
            max_frames=metric_max_frames,
        )
        if not struct_df.empty:
            struct_df.to_csv(output_dir / "structural_metrics.csv", index=False)
            ddg_df = merge_with_structural_metrics(ddg_df, struct_df)

    logging.info("Computing boundness QC")
    qc_df = compute_boundness_qc(run_df, ligand_resname)
    if not qc_df.empty:
        qc_df.to_csv(output_dir / "boundness_qc.csv", index=False)

    logging.info("Computing RMSD convergence profiles")
    rmsd_df = collect_ca_rmsd_profiles(run_df)
    if not rmsd_df.empty:
        rmsd_df.to_csv(output_dir / "rmsd_ca_profiles.csv", index=False)

    logging.info("Computing DOR-RT COM distance convergence profiles")
    com_df = collect_com_distance_profiles(
        run_df,
        ligand_resname=ligand_resname,
        frame_stride=metric_frame_stride,
        max_frames=max(400, metric_max_frames),
    )
    if not com_df.empty:
        com_df.to_csv(output_dir / "com_distance_profiles.csv", index=False)

    ddg_df.to_csv(output_dir / "ddg_full.csv", index=False)

    correlations = compute_correlations(ddg_df)
    correlations.to_csv(output_dir / "correlation_analysis.csv", index=False)

    return correlations, ddg_df
