#!/usr/bin/env python3
"""Recompute 100 ns MD convergence descriptors from analysis DCDs.

Does not reuse cached RMSD/COM tables. Loads pre-imaged *_analysis_pbcfix.dcd
(write those first with fix_pbc_trajectories). Library calls:

  - protein Cα RMSD vs frame 0          (mdtraj.rmsd)
  - DOR heavy-atom RMSD vs frame 0      (mdtraj.rmsd after NNIBP Cα fit)
  - DOR–RT COM distance                 (min-image, after PBC)
  - reporter min heavy-atom distances   (numpy; Ser105, Val179, residue 227)

Convergence figures are raw coordinate-vs-time traces (faint per-rep + mean),
the same style as the old RMSD / COM panels. Not last-N window tests.

Time axis comes from infer_production_ns() (DCD fingerprint / repaired JSON), not DCD dt.

Example:
    ~/miniconda3/envs/nnrti-prep/bin/python -m src.analysis.cli.compute_md_convergence
"""
from __future__ import annotations

import argparse
import json
import logging
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ..md_timing import infer_production_ns
from ..pbc import load_mdtraj_trajectory, pbcfix_dcd_for, raw_analysis_dcd_for
from ..result_collector import _prepare_profile_jobs, collect_md_results

REPO = Path(__file__).resolve().parents[3]
LOGGER = logging.getLogger("md_convergence")

# Canonical HIV-1 RT numbering. traj resSeq = auth + resid_offset.
RESID_OFFSET = -3
REPORTERS = {
    "ser105": 105,
    "val179": 179,
    "res227": 227,
}
NNIBP_AUTH = (100, 101, 103, 106, 107, 108, 179, 181, 188, 189, 190, 227, 229, 234, 318)
CONTACT_CUTOFF_A = 4.0
METRIC_LABELS = {
    "ca_rmsd_angstrom": "Cα RMSD (Å)",
    "dor_rmsd_angstrom": "DOR pose RMSD (Å)",
    "com_distance_angstrom": "DOR–RT COM distance (Å)",
    "ser105_min_angstrom": "Ser105–DOR min distance (Å)",
    "val179_min_angstrom": "Val179–DOR min distance (Å)",
    "res227_min_angstrom": "Residue 227–DOR min distance (Å)",
    "nnibp_rg_angstrom": "NNIBP Cα Rg (Å)",
}
TRACE_COLORS = {
    "ca_rmsd_angstrom": "#1f77b4",
    "dor_rmsd_angstrom": "#2ca02c",
    "com_distance_angstrom": "#d62728",
    "ser105_min_angstrom": "#9467bd",
    "val179_min_angstrom": "#8c564b",
    "res227_min_angstrom": "#e377c2",
    "nnibp_rg_angstrom": "#7f7f7f",
}


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Recompute MD convergence traces from analysis DCDs.")
    parser.add_argument("--manifest", type=Path, default=Path("manifests/md_manifest.csv"))
    parser.add_argument("--output-dir", type=Path, default=Path("results/analysis/md_convergence"))
    parser.add_argument("--resid-offset", type=int, default=RESID_OFFSET)
    parser.add_argument("--contact-cutoff", type=float, default=CONTACT_CUTOFF_A)
    parser.add_argument("--n-blocks", type=int, default=5)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--mutations", nargs="*", default=None, help="Subset of mutation labels (e.g. WT V106A)")
    parser.add_argument(
        "--from-tables",
        action="store_true",
        help="Replot from existing frame_traces.csv; do not reload trajectories.",
    )
    return parser.parse_args()


def _display_mutation(raw: str) -> str:
    text = str(raw).strip()
    return "WT" if text.lower() == "wt" else text


def _select_resseq(topology, resseq: int, *, heavy: bool) -> np.ndarray:
    query = f"chainid 0 and resSeq {int(resseq)}"
    if heavy:
        query += " and not element H"
    idx = np.asarray(topology.select(query), dtype=int)
    return idx


def _min_heavy_distance_nm(xyz: np.ndarray, idx_a: np.ndarray, idx_b: np.ndarray) -> np.ndarray:
    a = xyz[:, idx_a, :]
    b = xyz[:, idx_b, :]
    delta = a[:, :, None, :] - b[:, None, :, :]
    return np.sqrt(np.sum(delta * delta, axis=-1)).min(axis=(1, 2))


def _process_job(job: dict, resid_offset: int, contact_cutoff: float) -> tuple[pd.DataFrame, dict]:
    import mdtraj as md

    mutation = _display_mutation(job["mutation"])
    replicate = int(job["replicate"])
    production_ps = float(job["production_ps"])
    traj = load_mdtraj_trajectory(Path(job["trajectory"]), Path(job["topology"]))

    ca_idx = np.asarray(traj.topology.select("protein and name CA"), dtype=int)
    lig_idx = np.asarray(traj.topology.select("resname '2KW' and not element H"), dtype=int)
    prot_idx = np.asarray(traj.topology.select("protein"), dtype=int)
    lig_all = np.asarray(traj.topology.select("resname '2KW'"), dtype=int)
    if ca_idx.size == 0:
        raise ValueError("no protein Cα atoms")
    if lig_idx.size == 0:
        raise ValueError("no DOR heavy atoms")
    if prot_idx.size == 0 or lig_all.size == 0:
        raise ValueError("empty protein or ligand selection for COM")

    nnibp_resseq = [int(auth) + int(resid_offset) for auth in NNIBP_AUTH]
    nnibp_query = " or ".join(f"resSeq {r}" for r in nnibp_resseq)
    nnibp_ca = np.asarray(
        traj.topology.select(f"chainid 0 and name CA and ({nnibp_query})"),
        dtype=int,
    )
    if nnibp_ca.size < 8:
        raise ValueError(f"NNIBP Cα selection too small ({nnibp_ca.size})")

    reporter_idx = {}
    for key, auth in REPORTERS.items():
        idx = _select_resseq(traj.topology, int(auth) + int(resid_offset), heavy=True)
        if idx.size == 0:
            raise ValueError(f"no heavy atoms for {key} (auth {auth})")
        reporter_idx[key] = idx

    n_frames = int(traj.n_frames)
    if n_frames < 2:
        raise ValueError("need at least 2 frames")
    time_ns = np.arange(n_frames, dtype=float) * (production_ps / (n_frames - 1)) / 1000.0

    prot_com = traj.xyz[:, prot_idx, :].mean(axis=1)
    lig_com = traj.xyz[:, lig_all, :].mean(axis=1)
    delta = lig_com - prot_com
    if traj.unitcell_lengths is not None:
        box = np.asarray(traj.unitcell_lengths, dtype=float)
        delta = delta - box * np.round(delta / box)
    com_dist = np.linalg.norm(delta, axis=1) * 10.0

    traj.superpose(traj, frame=0, atom_indices=ca_idx)
    ca_rmsd = md.rmsd(traj, traj, frame=0, atom_indices=ca_idx) * 10.0

    pose = traj[:]
    pose.superpose(pose, frame=0, atom_indices=nnibp_ca)
    dor_rmsd = md.rmsd(pose, pose, frame=0, atom_indices=lig_idx) * 10.0
    nnibp_rg = md.compute_rg(traj.atom_slice(nnibp_ca)) * 10.0

    rows = []
    distances = {}
    for key, idx in reporter_idx.items():
        dist_a = _min_heavy_distance_nm(traj.xyz, idx, lig_idx) * 10.0
        distances[key] = dist_a
    for i in range(n_frames):
        row = {
            "mutation": mutation,
            "safe_label": str(job["safe_label"]),
            "replicate": replicate,
            "frame_index": i,
            "time_ns": float(time_ns[i]),
            "total_ns": float(production_ps / 1000.0),
            "n_frames": n_frames,
            "ca_rmsd_angstrom": float(ca_rmsd[i]),
            "dor_rmsd_angstrom": float(dor_rmsd[i]),
            "com_distance_angstrom": float(com_dist[i]),
            "nnibp_rg_angstrom": float(nnibp_rg[i]),
        }
        for key, dist in distances.items():
            row[f"{key}_min_angstrom"] = float(dist[i])
            row[f"{key}_occupancy"] = float(dist[i] < contact_cutoff)
        rows.append(row)

    meta = {
        "mutation": mutation,
        "replicate": replicate,
        "n_frames": n_frames,
        "total_ns": float(production_ps / 1000.0),
        "trajectory": str(job["trajectory"]),
        "n_ca": int(ca_idx.size),
        "n_dor_heavy": int(lig_idx.size),
        "n_nnibp_ca": int(nnibp_ca.size),
        "timing_source": str(job.get("timing_source") or ""),
    }
    return pd.DataFrame(rows), meta


def _running_mean(values: np.ndarray) -> np.ndarray:
    return np.cumsum(values, dtype=float) / np.arange(1, values.size + 1)


def _interp_on_grid(time_ns: np.ndarray, values: np.ndarray, grid: np.ndarray) -> np.ndarray:
    finite = np.isfinite(time_ns) & np.isfinite(values)
    if finite.sum() < 2:
        return np.full(grid.shape, np.nan)
    order = np.argsort(time_ns[finite])
    t = time_ns[finite][order]
    y = values[finite][order]
    t_unique, idx = np.unique(t, return_index=True)
    return np.interp(grid, t_unique, y[idx], left=np.nan, right=np.nan)


def _mean_sem(values: np.ndarray) -> tuple[float, float]:
    vals = values[np.isfinite(values)]
    if vals.size == 0:
        return float("nan"), float("nan")
    mean = float(vals.mean())
    sem = float(vals.std(ddof=1) / np.sqrt(vals.size)) if vals.size > 1 else 0.0
    return mean, sem


def _window_stats(frame_df: pd.DataFrame, metrics: list[str]) -> pd.DataFrame:
    rows = []
    for (mutation, replicate), sub in frame_df.groupby(["mutation", "replicate"], sort=False):
        sub = sub.sort_values("time_ns")
        t = sub["time_ns"].to_numpy(dtype=float)
        tmax = float(np.nanmax(t)) if t.size else float("nan")
        windows = {
            "full": np.ones(t.size, dtype=bool),
            "first_20ns": t <= 20.0,
            "last_20ns": t >= (tmax - 20.0),
            "last_50ns": t >= (tmax - 50.0),
        }
        for window, mask in windows.items():
            if mask.sum() < 2:
                continue
            row = {
                "mutation": mutation,
                "replicate": int(replicate),
                "window": window,
                "n_frames": int(mask.sum()),
                "t_start_ns": float(t[mask].min()),
                "t_end_ns": float(t[mask].max()),
                "total_ns": float(sub["total_ns"].iloc[0]),
            }
            for metric in metrics:
                mean, sem = _mean_sem(sub.loc[mask, metric].to_numpy(dtype=float))
                row[f"{metric}_mean"] = mean
                row[f"{metric}_sem"] = sem
            rows.append(row)
    return pd.DataFrame(rows)


def _block_stats(frame_df: pd.DataFrame, metrics: list[str], n_blocks: int) -> pd.DataFrame:
    rows = []
    for (mutation, replicate), sub in frame_df.groupby(["mutation", "replicate"], sort=False):
        sub = sub.sort_values("time_ns")
        t = sub["time_ns"].to_numpy(dtype=float)
        if t.size < n_blocks * 2:
            continue
        edges = np.linspace(float(t.min()), float(t.max()), n_blocks + 1)
        for block_id in range(n_blocks):
            lo, hi = edges[block_id], edges[block_id + 1]
            mask = (t >= lo) & (t <= hi if block_id == n_blocks - 1 else t < hi)
            if mask.sum() < 2:
                continue
            row = {
                "mutation": mutation,
                "replicate": int(replicate),
                "block_id": block_id + 1,
                "n_blocks": n_blocks,
                "t_start_ns": float(t[mask].min()),
                "t_end_ns": float(t[mask].max()),
                "n_frames": int(mask.sum()),
            }
            for metric in metrics:
                mean, sem = _mean_sem(sub.loc[mask, metric].to_numpy(dtype=float))
                row[f"{metric}_mean"] = mean
                row[f"{metric}_sem"] = sem
            rows.append(row)
    return pd.DataFrame(rows)


def _genotype_running_mean(frame_df: pd.DataFrame, metrics: list[str], n_grid: int = 200) -> pd.DataFrame:
    rows = []
    for mutation, mut_df in frame_df.groupby("mutation", sort=False):
        tmax = float(mut_df["time_ns"].max())
        grid = np.linspace(0.0, tmax, n_grid)
        for metric in metrics:
            stacked = []
            for _, sub in mut_df.groupby("replicate"):
                sub = sub.sort_values("time_ns")
                run = _running_mean(sub[metric].to_numpy(dtype=float))
                stacked.append(_interp_on_grid(sub["time_ns"].to_numpy(dtype=float), run, grid))
            arr = np.vstack(stacked)
            mean = np.nanmean(arr, axis=0)
            n = np.sum(np.isfinite(arr), axis=0)
            std = np.nanstd(arr, axis=0, ddof=1)
            sem = np.where(n > 1, std / np.sqrt(n), 0.0)
            for i, t in enumerate(grid):
                rows.append(
                    {
                        "mutation": mutation,
                        "metric": metric,
                        "time_ns": float(t),
                        "running_mean": float(mean[i]),
                        "running_sem": float(sem[i]),
                        "n_replicates": int(n[i]),
                    }
                )
    return pd.DataFrame(rows)


def _mutation_sort_key(mutation: str) -> tuple[int, str]:
    if mutation == "WT":
        return (0, mutation)
    if "+" in mutation:
        return (2, mutation)
    return (1, mutation)


def _interp_mean_trace(sub: pd.DataFrame, y_col: str, n_grid: int = 200) -> tuple[np.ndarray | None, np.ndarray | None]:
    reps = []
    for _, rep_df in sub.groupby("replicate"):
        rep_df = rep_df.sort_values("time_ns")
        x = rep_df["time_ns"].to_numpy(dtype=float)
        y = rep_df[y_col].to_numpy(dtype=float)
        keep = np.isfinite(x) & np.isfinite(y)
        x, y = x[keep], y[keep]
        if x.size < 2:
            continue
        order = np.argsort(x)
        x, y = x[order], y[order]
        keep = np.r_[True, np.diff(x) > 0]
        x, y = x[keep], y[keep]
        if x.size < 2:
            continue
        reps.append((x, y))
    if not reps:
        return None, None
    grid = np.linspace(min(x.min() for x, _ in reps), max(x.max() for x, _ in reps), n_grid)
    stacked = []
    for x, y in reps:
        yi = np.interp(grid, x, y)
        yi[(grid < x.min()) | (grid > x.max())] = np.nan
        stacked.append(yi)
    return grid, np.nanmean(np.vstack(stacked), axis=0)


def _plot_coordinate_traces(
    frame_df: pd.DataFrame,
    y_col: str,
    out_path: Path,
    *,
    title: str,
    ylabel: str,
    color: str,
) -> None:
    if y_col not in frame_df.columns:
        LOGGER.warning("skip %s: column %s missing", out_path.name, y_col)
        return
    muts = sorted(frame_df["mutation"].dropna().unique(), key=_mutation_sort_key)
    if not muts:
        return
    max_cols = 4
    ncols = min(len(muts), max_cols)
    nrows = int(np.ceil(len(muts) / max_cols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(5.5 * ncols, 3.5 * nrows), squeeze=False)
    for plot_i, mutation in enumerate(muts):
        ax = axes[plot_i // max_cols][plot_i % max_cols]
        sub = frame_df[frame_df["mutation"] == mutation]
        for _, rep_df in sub.groupby("replicate"):
            rep_df = rep_df.sort_values("time_ns")
            ax.plot(rep_df["time_ns"], rep_df[y_col], linewidth=1.0, alpha=0.5, color=color)
        x_mean, y_mean = _interp_mean_trace(sub, y_col)
        if x_mean is not None:
            ax.plot(x_mean, y_mean, linewidth=2.0, color=color, label="replicate mean")
        ax.set_title(str(mutation), fontsize=10, fontweight="bold")
        if plot_i % max_cols == 0:
            ax.set_ylabel(ylabel, fontsize=9)
        ax.set_xlabel("Time (ns)", fontsize=9)
        ax.grid(alpha=0.2, linestyle=":")
    for ax in axes.ravel()[len(muts) :]:
        ax.set_visible(False)
    fig.suptitle(title, y=0.995, fontsize=14, fontweight="bold")
    fig.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def _write_convergence_plots(frame_df: pd.DataFrame, plots: Path) -> None:
    for stale in plots.glob("*.png"):
        stale.unlink()
    panels = [
        (
            "ca_rmsd_angstrom",
            "rmsd_convergence.png",
            "Cα RMSD from the first production frame",
        ),
        (
            "dor_rmsd_angstrom",
            "dor_rmsd_convergence.png",
            "DOR pose RMSD after NNIBP Cα fit",
        ),
        (
            "com_distance_angstrom",
            "com_distance_convergence.png",
            "DOR–RT center-of-mass distance",
        ),
        (
            "ser105_min_angstrom",
            "ser105_distance_convergence.png",
            "Ser105–DOR minimum heavy-atom distance",
        ),
        (
            "val179_min_angstrom",
            "val179_distance_convergence.png",
            "Val179–DOR minimum heavy-atom distance",
        ),
        (
            "res227_min_angstrom",
            "res227_distance_convergence.png",
            "Residue 227–DOR minimum heavy-atom distance",
        ),
    ]
    for y_col, filename, title in panels:
        _plot_coordinate_traces(
            frame_df,
            y_col,
            plots / filename,
            title=title,
            ylabel=METRIC_LABELS[y_col],
            color=TRACE_COLORS[y_col],
        )


def main() -> int:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    args = _parse_args()
    out = args.output_dir
    tables = out / "tables"
    plots = out / "plots"
    config = out / "config"
    for folder in (tables, plots, config):
        folder.mkdir(parents=True, exist_ok=True)

    if args.from_tables:
        frame_path = tables / "frame_traces.csv"
        if not frame_path.exists():
            raise SystemExit(f"--from-tables requires {frame_path}")
        frame_df = pd.read_csv(frame_path)
        _write_convergence_plots(frame_df, plots)
        LOGGER.info("Replotted %d frames from %s", len(frame_df), frame_path)
        return 0

    run_df = collect_md_results(args.manifest)
    if args.mutations:
        wanted = {_display_mutation(m) for m in args.mutations}
        run_df = run_df[run_df["mutation"].map(_display_mutation).isin(wanted)].copy()
    jobs = _prepare_profile_jobs(run_df)
    if not jobs:
        raise SystemExit("No analysis trajectories found.")

    missing_pbc = []
    for job in jobs:
        raw = Path(job["trajectory"])
        imaged = pbcfix_dcd_for(raw) if not raw.name.endswith("_pbcfix.dcd") else raw
        if not imaged.exists():
            missing_pbc.append(str(imaged))
            continue
        job["trajectory"] = str(imaged)
    if missing_pbc:
        preview = "\n  ".join(missing_pbc[:8])
        raise SystemExit(
            f"Missing {len(missing_pbc)} PBC-imaged DCDs. Write them first:\n"
            f"  python -m src.analysis.cli.fix_pbc_trajectories --workers 4\n"
            f"  {preview}"
        )

    timing_rows = []
    for job in jobs:
        traj_path = Path(job["trajectory"])
        raw_dcd = raw_analysis_dcd_for(traj_path)
        timing_dcd = raw_dcd if raw_dcd.exists() else traj_path
        safe = str(job["safe_label"])
        rep = int(job["replicate"])
        call = infer_production_ns(
            dcd_path=timing_dcd,
            json_path=traj_path.parent / f"{safe}_rep{rep:02d}.json",
            state_csv_path=traj_path.parent / f"{safe}_rep{rep:02d}_md_state.csv",
            mutation=_display_mutation(job["mutation"]),
            replicate=rep,
        )
        job["production_ps"] = float(call.production_ns) * 1000.0
        job["timing_source"] = call.source
        timing_rows.append(call.__dict__)
        if call.note:
            LOGGER.info(
                "TIME %s rep%s: %.1f ns via %s (%s)",
                call.mutation,
                call.replicate,
                call.production_ns,
                call.source,
                call.note,
            )
    pd.DataFrame(timing_rows).to_csv(tables / "timing_audit.csv", index=False)

    metrics = [
        "ca_rmsd_angstrom",
        "dor_rmsd_angstrom",
        "com_distance_angstrom",
        "ser105_min_angstrom",
        "val179_min_angstrom",
        "res227_min_angstrom",
        "nnibp_rg_angstrom",
        "ser105_occupancy",
        "val179_occupancy",
        "res227_occupancy",
    ]

    frames: list[pd.DataFrame] = []
    metas: list[dict] = []
    failures: list[dict] = []
    LOGGER.info("Processing %d trajectories with %d workers", len(jobs), args.workers)
    with ProcessPoolExecutor(max_workers=max(1, int(args.workers))) as pool:
        futs = {
            pool.submit(_process_job, job, int(args.resid_offset), float(args.contact_cutoff)): job
            for job in jobs
        }
        done = 0
        for fut in as_completed(futs):
            job = futs[fut]
            done += 1
            try:
                frame_df, meta = fut.result()
                frames.append(frame_df)
                metas.append(meta)
                LOGGER.info(
                    "[%d/%d] %s rep%s  %d frames  %.1f ns",
                    done,
                    len(jobs),
                    meta["mutation"],
                    meta["replicate"],
                    meta["n_frames"],
                    meta["total_ns"],
                )
            except Exception as exc:
                failures.append(
                    {
                        "mutation": _display_mutation(job["mutation"]),
                        "replicate": int(job["replicate"]),
                        "trajectory": str(job["trajectory"]),
                        "error": str(exc),
                    }
                )
                LOGGER.warning(
                    "[%d/%d] FAILED %s rep%s: %s",
                    done,
                    len(jobs),
                    job["mutation"],
                    job["replicate"],
                    exc,
                )

    if not frames:
        raise SystemExit("All trajectories failed.")

    frame_df = pd.concat(frames, ignore_index=True)
    meta_df = pd.DataFrame(metas)
    fail_df = pd.DataFrame(failures)
    run_mean_df = _genotype_running_mean(frame_df, metrics)
    window_df = _window_stats(frame_df, metrics)
    block_df = _block_stats(frame_df, metrics, int(args.n_blocks))

    frame_df.to_csv(tables / "frame_traces.csv", index=False)
    meta_df.to_csv(tables / "replicate_inventory.csv", index=False)
    run_mean_df.to_csv(tables / "running_means.csv", index=False)
    window_df.to_csv(tables / "window_comparison.csv", index=False)
    block_df.to_csv(tables / "block_averages.csv", index=False)
    if not fail_df.empty:
        fail_df.to_csv(tables / "failures.csv", index=False)

    _write_convergence_plots(frame_df, plots)

    config_payload = {
        "manifest": str(args.manifest),
        "resid_offset": int(args.resid_offset),
        "contact_cutoff_angstrom": float(args.contact_cutoff),
        "n_blocks": int(args.n_blocks),
        "reporters_auth": REPORTERS,
        "nnibp_auth": list(NNIBP_AUTH),
        "n_jobs": len(jobs),
        "n_success": int(len(metas)),
        "n_failed": int(len(failures)),
        "libraries": ["mdtraj", "numpy", "matplotlib"],
        "notes": [
            "Convergence figures are raw coordinate-vs-time traces (per-rep + mean), not last-N windows.",
            "Loads pre-imaged *_analysis_pbcfix.dcd; does not PBC-correct on the fly.",
            "RMSD is vs frame 0, not a cached table.",
            "DOR pose RMSD uses NNIBP Cα superposition then ligand heavy-atom RMSD.",
            "COM is protein-all vs ligand-all, min-image.",
            "Time axis uses infer_production_ns() on the raw analysis DCD.",
        ],
    }
    (config / "run_config.json").write_text(json.dumps(config_payload, indent=2) + "\n")
    LOGGER.info("Wrote %s (%d frames, %d successes, %d failures)", out, len(frame_df), len(metas), len(failures))
    return 0 if not failures else 1


if __name__ == "__main__":
    raise SystemExit(main())
