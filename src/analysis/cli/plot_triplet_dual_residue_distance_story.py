#!/usr/bin/env python3
"""Plot a 100 ns two-panel residue-to-DOR distance triplet story on a common time grid."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .plot_triplet_residue_distance_story import (
    MUTATION_COLORS,
    _build_common_time_grid,
    _compute_residue_trace,
    _interpolate_replicates_to_grid,
    _load_fold_map,
    _load_replicate_meta,
)


def _parse_csv_tokens(text: str) -> list[str]:
    return [token.strip() for token in str(text).split(",") if token.strip()]


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plot a two-panel residue-to-DOR triplet story over 100 ns.")
    parser.add_argument(
        "--manifest",
        type=Path,
        default=Path("manifests/md_manifest.csv"),
    )
    parser.add_argument(
        "--susceptibility-xlsx",
        type=Path,
        default=Path("data/DRM-susceptibilities.csv.xlsx"),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results/analysis/triplet_story_analyses/geometry_wt_v106a_f227l_v106i_f227c_two_panel"),
    )
    parser.add_argument("--triplet", type=str, default="WT,V106A+F227L,V106I+F227C")
    parser.add_argument("--auth-resseqs", type=str, default="105,227")
    parser.add_argument("--metric-names", type=str, default="SER105-DOR,Pos227-DOR")
    parser.add_argument(
        "--ylabels",
        type=str,
        default="Min SER105-DOR Distance (A),Min residue 227-DOR Distance (A)",
    )
    parser.add_argument("--titles", type=str, default="SER105,Position 227")
    parser.add_argument("--output-prefix", type=str, default="triplet_story_100ns_WT_V106A_F227L_V106I_F227C_SER105_POS227")
    parser.add_argument("--max-time-ns", type=float, default=100.0)
    parser.add_argument(
        "--force-total-ns",
        type=float,
        default=None,
        help="If set, ignore per-replicate timing metadata and map each trajectory uniformly onto this duration.",
    )
    parser.add_argument("--resid-offset", type=int, default=-3)
    parser.add_argument("--ligand-resname", type=str, default="2KW")
    parser.add_argument(
        "--triplet-colors",
        type=str,
        default="",
        help="Optional comma-separated colors to apply to the triplet in order.",
    )
    return parser.parse_args()


def main() -> int:
    args = _parse_args()
    if not args.manifest.exists():
        raise FileNotFoundError(args.manifest)
    if not args.susceptibility_xlsx.exists():
        raise FileNotFoundError(args.susceptibility_xlsx)

    triplet = _parse_csv_tokens(args.triplet)
    if len(triplet) != 3:
        raise ValueError("--triplet must contain exactly three comma-separated mutations")

    auth_resseqs = [int(token) for token in _parse_csv_tokens(args.auth_resseqs)]
    metric_names = _parse_csv_tokens(args.metric_names)
    ylabels = _parse_csv_tokens(args.ylabels)
    titles = _parse_csv_tokens(args.titles)
    n_panels = len(auth_resseqs)
    if n_panels < 1:
        raise ValueError("At least one auth residue is required")
    if not (len(metric_names) == len(ylabels) == len(titles) == n_panels):
        raise ValueError("--auth-resseqs, --metric-names, --ylabels, and --titles must have equal lengths")

    triplet_colors = _parse_csv_tokens(args.triplet_colors)
    if triplet_colors and len(triplet_colors) != len(triplet):
        raise ValueError("--triplet-colors must provide exactly one color per triplet entry")

    out_tables = args.output_dir / "tables"
    out_plots = args.output_dir / "plots"
    out_config = args.output_dir / "config"
    out_tables.mkdir(parents=True, exist_ok=True)
    out_plots.mkdir(parents=True, exist_ok=True)
    out_config.mkdir(parents=True, exist_ok=True)

    fold_map = _load_fold_map(args.susceptibility_xlsx)
    metas = _load_replicate_meta(args.manifest, needed_mutations=set(triplet))
    if not metas:
        raise ValueError("No valid replicate metadata found for requested triplet")

    raw_metric_frames: list[pd.DataFrame] = []
    interp_metric_frames: list[pd.DataFrame] = []
    mean_rows: list[pd.DataFrame] = []

    for auth_resseq, metric_name in zip(auth_resseqs, metric_names):
        metric_rows: list[pd.DataFrame] = []
        for meta in sorted(metas, key=lambda x: (x.mutation, x.replicate)):
            t_sel, trace = _compute_residue_trace(
                analysis_dcd=str(meta.analysis_dcd),
                topology_pdb=str(meta.topology_pdb),
                total_ns=float(args.force_total_ns) if args.force_total_ns is not None else float(meta.total_ns),
                window_ns=float(args.max_time_ns),
                auth_resseq=int(auth_resseq),
                resid_offset=int(args.resid_offset),
                ligand_resname=str(args.ligand_resname),
            )
            sub = pd.DataFrame(
                {
                    "mutation": str(meta.mutation),
                    "replicate": int(meta.replicate),
                    "time_ns": t_sel.astype(float),
                    "metric": str(metric_name),
                    "metric_value": trace,
                }
            )
            metric_rows.append(sub)
        metric_df = pd.concat(metric_rows, ignore_index=True)
        metric_df = metric_df.sort_values(["mutation", "replicate", "time_ns"], kind="stable").reset_index(drop=True)
        raw_metric_frames.append(metric_df)

        time_grid = _build_common_time_grid(metric_df, max_time_ns=float(args.max_time_ns))
        interp_df = _interpolate_replicates_to_grid(
            metric_df.rename(columns={"metric_value": metric_name})[
                ["mutation", "replicate", "time_ns", metric_name]
            ],
            metric_name,
            time_grid,
        )
        interp_df["metric"] = str(metric_name)
        interp_metric_frames.append(interp_df.rename(columns={metric_name: "metric_value"}))

        mean_df = (
            interp_df.groupby(["mutation", "time_ns"], as_index=False)[metric_name]
            .agg(["mean", "std", "count"])
            .reset_index()
            .rename(columns={"mean": "metric_mean", "std": "metric_std", "count": "n_replicates"})
        )
        mean_df["metric"] = str(metric_name)
        mean_df["metric_sem"] = (
            mean_df["metric_std"].fillna(0.0) / mean_df["n_replicates"].clip(lower=1).pow(0.5)
        ).astype(float)
        mean_rows.append(
            mean_df[["mutation", "time_ns", "metric", "metric_mean", "metric_std", "metric_sem", "n_replicates"]]
        )

    raw_df_all = pd.concat(raw_metric_frames, ignore_index=True)
    interp_df_all = pd.concat(interp_metric_frames, ignore_index=True)
    mean_df_all = pd.concat(mean_rows, ignore_index=True)
    raw_df_all.to_csv(out_tables / "trace_values_raw_metrics.csv", index=False)
    interp_df_all.to_csv(out_tables / "trace_values.csv", index=False)
    mean_df_all.to_csv(out_tables / "mean_traces.csv", index=False)

    fig, axes = plt.subplots(1, n_panels, figsize=(7.0 * n_panels, 5.6), sharex=True, constrained_layout=True)
    if n_panels == 1:
        axes = [axes]

    legend_handles = None
    legend_labels = None
    xmax = 0.0
    for ax, metric_name, ylabel, title in zip(axes, metric_names, ylabels, titles):
        sub = mean_df_all[mean_df_all["metric"].astype(str) == str(metric_name)].copy()
        for idx, mutation in enumerate(triplet):
            color = triplet_colors[idx] if triplet_colors else MUTATION_COLORS.get(mutation, "#555555")
            mut_mean = sub[sub["mutation"].astype(str) == mutation].copy()
            x = mut_mean["time_ns"].to_numpy(dtype=float)
            y = mut_mean["metric_mean"].to_numpy(dtype=float)
            sem = mut_mean["metric_sem"].to_numpy(dtype=float)
            xmax = max(xmax, float(np.nanmax(x)) if len(x) else 0.0)
            fold = fold_map.get(mutation, float("nan"))
            label = f"{mutation} ({fold:.1f}x)" if pd.notna(fold) else str(mutation)
            ax.plot(x, y, color=color, linewidth=2.1, alpha=0.95, label=label)
            lo = y - sem
            hi = y + sem
            ok = np.isfinite(lo) & np.isfinite(hi)
            if np.any(ok):
                ax.fill_between(x[ok], lo[ok], hi[ok], color=color, alpha=0.16, linewidth=0)
        ax.axhline(4.0, color="#666666", linestyle=":", linewidth=1.1)
        ax.set_xlim(0.0, xmax if xmax > 0 else float(args.max_time_ns))
        ax.set_xlabel("Time (ns)")
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.grid(alpha=0.22, linestyle=":")
        for spine in ("top", "right"):
            ax.spines[spine].set_visible(False)
        if legend_handles is None:
            legend_handles, legend_labels = ax.get_legend_handles_labels()

    if legend_handles and legend_labels:
        legend_handles = list(legend_handles) + [plt.Line2D([0], [0], color="#666666", linestyle=":", linewidth=1.1)]
        legend_labels = list(legend_labels) + ["4.0 A reference"]
        axes[-1].legend(legend_handles, legend_labels, loc="upper right", frameon=True, fontsize=9)

    png = out_plots / f"{str(args.output_prefix)}.png"
    fig.savefig(png, dpi=300, bbox_inches="tight")
    plt.close(fig)

    (out_config / "run_config.json").write_text(
        json.dumps(
            {
                "manifest": str(args.manifest),
                "susceptibility_xlsx": str(args.susceptibility_xlsx),
                "output_dir": str(args.output_dir),
                "triplet": triplet,
                "max_time_ns": float(args.max_time_ns),
                "resid_offset": int(args.resid_offset),
                "ligand_resname": str(args.ligand_resname),
                "auth_resseqs": auth_resseqs,
                "metric_names": metric_names,
                "ylabels": ylabels,
                "titles": titles,
                "output_prefix": str(args.output_prefix),
                "triplet_colors": triplet_colors,
                "force_total_ns": None if args.force_total_ns is None else float(args.force_total_ns),
                "aggregation": "replicates interpolated onto common time grid before mean/SEM",
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
