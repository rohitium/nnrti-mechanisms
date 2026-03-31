#!/usr/bin/env python3
"""Plot a 100 ns triplet story for the PRO225-to-DOR distance."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ..susceptibility import load_dor_susceptibilities


MUTATION_COLORS = {
    "WT": "#333333",
    "Y181C": "#4c78a8",
    "Y188L": "#e45756",
}

METRIC = "residue_min_distance_PRO225_angstrom"
YLABEL = "Min PRO225-DOR Distance (A)"
TITLE = "PRO225"


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plot a WT/comparator/DRM PRO225 distance story over 100 ns.")
    parser.add_argument(
        "--frame-feature-csv",
        type=Path,
        default=Path("results/analysis/ligand_pocket_features/tables/frame_features.csv"),
    )
    parser.add_argument(
        "--susceptibility-xlsx",
        type=Path,
        default=Path("data/DRM-susceptibilities.csv.xlsx"),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results/analysis/triplet_geometry_story_100ns"),
    )
    parser.add_argument("--triplet", type=str, default="WT,Y181C,Y188L")
    parser.add_argument("--max-time-ns", type=float, default=100.0)
    return parser.parse_args()


def _load_fold_map(xlsx_path: Path) -> dict[str, float]:
    df = load_dor_susceptibilities(xlsx_path)
    fold_map = {str(row["mutation"]): float(row["dor_fold_reduction"]) for _, row in df.iterrows()}
    fold_map["WT"] = 1.0
    return fold_map


def _build_common_time_grid(trace_df: pd.DataFrame, max_time_ns: float) -> np.ndarray:
    times = pd.to_numeric(trace_df["time_ns"], errors="coerce").dropna().to_numpy(dtype=float)
    times = np.unique(np.round(times, 6))
    times = times[(times >= 0.0) & (times <= float(max_time_ns))]
    if times.size == 0:
        return np.linspace(0.0, float(max_time_ns), num=101, dtype=float)
    if times[0] > 0.0:
        times = np.concatenate([np.array([0.0], dtype=float), times])
    if times[-1] < float(max_time_ns):
        times = np.concatenate([times, np.array([float(max_time_ns)], dtype=float)])
    return times


def _interpolate_replicates_to_grid(trace_df: pd.DataFrame, metric: str, time_grid: np.ndarray) -> pd.DataFrame:
    rows: list[pd.DataFrame] = []
    for (mutation, replicate), rep_df in trace_df.groupby(["mutation", "replicate"], sort=True):
        rep_df = rep_df.sort_values("time_ns", kind="stable").copy()
        x = pd.to_numeric(rep_df["time_ns"], errors="coerce").to_numpy(dtype=float)
        y = pd.to_numeric(rep_df[metric], errors="coerce").to_numpy(dtype=float)
        ok = np.isfinite(x) & np.isfinite(y)
        x = x[ok]
        y = y[ok]
        if x.size == 0:
            continue
        x_unique, unique_idx = np.unique(x, return_index=True)
        y_unique = y[unique_idx]
        y_interp = np.interp(time_grid, x_unique, y_unique)
        rows.append(
            pd.DataFrame(
                {
                    "mutation": str(mutation),
                    "replicate": int(replicate),
                    "time_ns": time_grid.astype(float),
                    metric: y_interp.astype(float),
                }
            )
        )
    if not rows:
        return pd.DataFrame(columns=["mutation", "replicate", "time_ns", metric])
    return pd.concat(rows, ignore_index=True)


def main() -> int:
    args = _parse_args()
    if not args.frame_feature_csv.exists():
        raise FileNotFoundError(args.frame_feature_csv)
    if not args.susceptibility_xlsx.exists():
        raise FileNotFoundError(args.susceptibility_xlsx)

    triplet = [token.strip() for token in str(args.triplet).split(",") if token.strip()]
    if len(triplet) != 3:
        raise ValueError("--triplet must contain exactly three comma-separated mutations")

    out_tables = args.output_dir / "tables"
    out_plots = args.output_dir / "plots"
    out_config = args.output_dir / "config"
    out_tables.mkdir(parents=True, exist_ok=True)
    out_plots.mkdir(parents=True, exist_ok=True)
    out_config.mkdir(parents=True, exist_ok=True)

    frame_df = pd.read_csv(args.frame_feature_csv)
    fold_map = _load_fold_map(args.susceptibility_xlsx)

    raw_trace_df = frame_df[
        frame_df["mutation"].astype(str).isin(triplet)
        & (pd.to_numeric(frame_df["time_ns"], errors="coerce") <= float(args.max_time_ns))
    ][["mutation", "replicate", "time_ns", METRIC]].copy()
    raw_trace_df["time_ns"] = pd.to_numeric(raw_trace_df["time_ns"], errors="coerce").astype(float)
    raw_trace_df = raw_trace_df.sort_values(["mutation", "replicate", "time_ns"], kind="stable").reset_index(drop=True)
    raw_trace_df.to_csv(out_tables / "trace_values_raw.csv", index=False)

    time_grid = _build_common_time_grid(raw_trace_df, max_time_ns=float(args.max_time_ns))
    trace_df = _interpolate_replicates_to_grid(raw_trace_df, METRIC, time_grid)
    trace_df.to_csv(out_tables / "trace_values.csv", index=False)

    mean_df = (
        trace_df.groupby(["mutation", "time_ns"], as_index=False)[METRIC]
        .agg(["mean", "std", "count"])
        .reset_index()
        .rename(columns={"mean": "metric_mean", "std": "metric_std", "count": "n_replicates"})
    )
    mean_df["metric_sem"] = (
        mean_df["metric_std"].fillna(0.0) / mean_df["n_replicates"].clip(lower=1).pow(0.5)
    ).astype(float)
    mean_df.to_csv(out_tables / "mean_traces.csv", index=False)

    fig, ax = plt.subplots(figsize=(13.8, 5.6), constrained_layout=True)
    xmax = 0.0
    for mutation in triplet:
        color = MUTATION_COLORS.get(mutation, "#555555")
        mut_mean = mean_df[mean_df["mutation"].astype(str) == mutation].copy()
        x = mut_mean["time_ns"].to_numpy(dtype=float)
        y = mut_mean["metric_mean"].to_numpy(dtype=float)
        sem = mut_mean["metric_sem"].to_numpy(dtype=float)
        xmax = max(xmax, float(np.nanmax(x)) if len(x) else 0.0)
        fold = fold_map.get(mutation, float("nan"))
        label = f"{mutation} ({fold:.1f}x)" if pd.notna(fold) else str(mutation)
        ax.plot(
            x,
            y,
            color=color,
            linewidth=2.1,
            alpha=0.95,
            label=label,
        )
        lo = y - sem
        hi = y + sem
        ok = np.isfinite(lo) & np.isfinite(hi)
        if np.any(ok):
            ax.fill_between(x[ok], lo[ok], hi[ok], color=color, alpha=0.16, linewidth=0)
    ax.axhline(4.0, color="#666666", linestyle=":", linewidth=1.1, label="4.0 A reference")
    ax.set_xlim(0.0, xmax if xmax > 0 else float(args.max_time_ns))
    ax.set_xlabel("Time (ns)")
    ax.set_ylabel(YLABEL)
    ax.set_title(TITLE)
    ax.grid(alpha=0.22, linestyle=":")
    ax.legend(loc="upper right", frameon=True, fontsize=9)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)

    png = out_plots / "triplet_story_100ns_WT_Y181C_Y188L_PRO225.png"
    fig.savefig(png, dpi=300, bbox_inches="tight")
    fig.savefig(png.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)

    (out_config / "run_config.json").write_text(
        json.dumps(
            {
                "frame_feature_csv": str(args.frame_feature_csv),
                "susceptibility_xlsx": str(args.susceptibility_xlsx),
                "output_dir": str(args.output_dir),
                "triplet": triplet,
                "max_time_ns": float(args.max_time_ns),
                "metric": METRIC,
                "aggregation": "replicates interpolated onto common time grid before mean/SEM",
                "common_time_grid_size": int(len(time_grid)),
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
