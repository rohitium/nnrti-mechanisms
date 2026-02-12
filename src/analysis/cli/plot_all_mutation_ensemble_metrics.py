#!/usr/bin/env python3
"""Plot structural metrics for all mutations with replicate aggregation.

Metrics:
1) Contacts
2) H-bonds
3) Binding pocket volume proxy

Aggregation rule:
1) Average over trajectory frames within each replicate (already in structural_metrics.csv).
2) Average replicate-level values within each mutation.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

_KJ_TO_KCAL = 1.0 / 4.184
_A3_TO_NM3 = 1.0 / 1000.0


def _mutation_sort_key(m: str) -> tuple[int, str]:
    if m == "WT":
        return (0, m)
    if "+" in m:
        return (2, m)
    return (1, m)


def _sem(values: pd.Series) -> float:
    vals = pd.to_numeric(values, errors="coerce").dropna()
    if len(vals) <= 1:
        return 0.0
    return float(vals.std(ddof=1) / np.sqrt(len(vals)))


def main() -> int:
    parser = argparse.ArgumentParser(description="Plot all-mutation ensemble metrics summary.")
    parser.add_argument("--structural-metrics", type=Path, default=Path("results/structural_metrics.csv"))
    parser.add_argument("--mmgbsa-metrics", type=Path, default=Path("results/mmgbsa_replicate_metrics.csv"))
    parser.add_argument(
        "--structural-metrics-checkpoint",
        type=Path,
        default=Path("results/.checkpoints/.checkpoint_structural_metrics.csv"),
    )
    parser.add_argument(
        "--mmgbsa-metrics-checkpoint",
        type=Path,
        default=Path("results/.checkpoints/.checkpoint_mmgbsa_replicate_metrics.csv"),
    )
    parser.add_argument("--output", type=Path, default=Path("results/plots/all_metrics.png"))
    parser.add_argument("--output-csv", type=Path, default=Path("results/ensemble_metrics_by_mutation.csv"))
    args = parser.parse_args()

    if not args.structural_metrics.exists() and not args.structural_metrics_checkpoint.exists():
        raise FileNotFoundError(
            f"Neither {args.structural_metrics} nor {args.structural_metrics_checkpoint} exists"
        )

    sm = pd.read_csv(args.structural_metrics) if args.structural_metrics.exists() else pd.DataFrame()
    if args.structural_metrics_checkpoint.exists():
        sm_ckpt = pd.read_csv(args.structural_metrics_checkpoint)
        if "error" in sm_ckpt.columns:
            sm_ckpt = sm_ckpt[sm_ckpt["error"].isna()].copy()
        overlap_keys = {"mutation", "replicate"}
        if sm.empty or (
            overlap_keys.issubset(sm.columns)
            and overlap_keys.issubset(sm_ckpt.columns)
            and sm_ckpt[["mutation", "replicate"]].drop_duplicates().shape[0]
            > sm[["mutation", "replicate"]].drop_duplicates().shape[0]
        ):
            sm = sm_ckpt.copy()

    required = {"mutation", "replicate", "contact_count", "hbond_count", "pocket_volume_proxy"}
    missing = sorted(required - set(sm.columns))
    if missing:
        raise ValueError(f"Structural metrics file is missing required columns: {missing}")

    mm = pd.read_csv(args.mmgbsa_metrics) if args.mmgbsa_metrics.exists() else pd.DataFrame()
    if args.mmgbsa_metrics_checkpoint.exists():
        mm_ckpt = pd.read_csv(args.mmgbsa_metrics_checkpoint)
        overlap_keys = {"mutation", "replicate"}
        if mm.empty or (
            overlap_keys.issubset(mm.columns)
            and overlap_keys.issubset(mm_ckpt.columns)
            and mm_ckpt[["mutation", "replicate"]].drop_duplicates().shape[0]
            > mm[["mutation", "replicate"]].drop_duplicates().shape[0]
        ):
            mm = mm_ckpt.copy()
    if "binding_dg" not in mm.columns:
        raise ValueError("MM/GBSA metrics file is missing required column: binding_dg")

    rep = sm[["mutation", "replicate", "contact_count", "hbond_count", "pocket_volume_proxy"]].copy()
    rep = rep.rename(
        columns={
            "contact_count": "contact_mean",
            "hbond_count": "hbond_mean",
            "pocket_volume_proxy": "pocket_volume_mean",
        }
    )
    dg = mm[["mutation", "replicate", "binding_dg"]].copy()
    dg["binding_dg_mean"] = pd.to_numeric(dg["binding_dg"], errors="coerce") * _KJ_TO_KCAL
    rep = rep.merge(dg[["mutation", "replicate", "binding_dg_mean"]], on=["mutation", "replicate"], how="inner")
    rep = rep.dropna(subset=["mutation", "replicate"])
    if rep.empty:
        raise ValueError("No valid overlapping replicate rows found in structural/MMGBSA metrics.")
    rep["pocket_volume_mean"] = pd.to_numeric(rep["pocket_volume_mean"], errors="coerce") * _A3_TO_NM3

    by_mut = (
        rep.groupby("mutation", as_index=False)
        .agg(
            n_replicates=("replicate", "nunique"),
            binding_dg_mean=("binding_dg_mean", "mean"),
            binding_dg_sem=("binding_dg_mean", _sem),
            contact_mean=("contact_mean", "mean"),
            contact_sem=("contact_mean", _sem),
            hbond_mean=("hbond_mean", "mean"),
            hbond_sem=("hbond_mean", _sem),
            pocket_volume_mean=("pocket_volume_mean", "mean"),
            pocket_volume_sem=("pocket_volume_mean", _sem),
        )
        .sort_values("mutation", key=lambda s: s.map(lambda m: _mutation_sort_key(str(m))))
        .reset_index(drop=True)
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    by_mut.to_csv(args.output_csv, index=False)

    order = by_mut["mutation"].tolist()
    x = np.arange(len(order), dtype=float)

    fig, axes = plt.subplots(4, 1, figsize=(max(12, 0.7 * len(order)), 13), sharex=True)
    metric_specs = [
        ("binding_dg_mean", "binding_dg_sem", "Binding Energy (kcal/mol)", "#9467bd"),
        ("contact_mean", "contact_sem", "RT-DOR Contacts (count)", "#1f77b4"),
        ("hbond_mean", "hbond_sem", "RT-DOR H-Bonds (count)", "#d62728"),
        ("pocket_volume_mean", "pocket_volume_sem", "Binding Pocket Volume (nm^3)", "#2ca02c"),
    ]

    rep_plot = rep.copy()
    rep_plot["mutation"] = pd.Categorical(rep_plot["mutation"], categories=order, ordered=True)
    rep_plot = rep_plot.sort_values(["mutation", "replicate"])

    for ax, (val_col, sem_col, ylabel, color) in zip(axes, metric_specs):
        metric_key = val_col.replace("_mean", "")
        # Light replicate points
        for i, mut in enumerate(order):
            sub = rep_plot[rep_plot["mutation"] == mut]
            if sub.empty:
                continue
            ycol = f"{metric_key}_mean"
            jitter = np.linspace(-0.12, 0.12, num=len(sub))
            ax.scatter(
                np.full(len(sub), i, dtype=float) + jitter,
                sub[ycol].to_numpy(dtype=float),
                color=color,
                alpha=0.5,
                s=26,
                linewidths=0.0,
                zorder=2,
            )

        ax.errorbar(
            x,
            by_mut[val_col].to_numpy(dtype=float),
            yerr=by_mut[sem_col].to_numpy(dtype=float),
            fmt="-o",
            color=color,
            linewidth=2.2,
            markersize=5,
            capsize=3,
            zorder=3,
            label="Mutation mean \u00b1 SEM across replicates",
        )
        wt_row = by_mut[by_mut["mutation"] == "WT"]
        if not wt_row.empty and pd.notna(wt_row[val_col].iloc[0]):
            wt_mean = float(wt_row[val_col].iloc[0])
            ax.axhline(
                wt_mean,
                color="#444444",
                linestyle="--",
                linewidth=1.0,
                alpha=0.9,
                label="WT mean",
            )
        ax.set_ylabel(ylabel)
        ax.grid(alpha=0.3, linestyle=":")
        # Legend intentionally omitted per figure style request.

    axes[0].set_title("Binding Metrics", fontsize=12, fontweight="bold")
    axes[-1].set_xticks(x)
    axes[-1].set_xticklabels(order, rotation=45, ha="right")
    axes[-1].set_xlabel("Mutation")

    fig.tight_layout()
    fig.savefig(args.output, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {args.output}")
    print(f"Wrote {args.output_csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
