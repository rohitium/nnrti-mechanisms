#!/usr/bin/env python3
"""Plot global mechanistic signatures vs doravirine fold-reduction.

Inputs:
  - results/analysis/binding_energy/tables/ddg_full.csv (replicate-level MM/GBSA + structural metrics)

Outputs:
  - results/plots/manuscript_global_signatures.png (default)
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def _sem(x: pd.Series) -> float:
    v = pd.to_numeric(x, errors="coerce").dropna().to_numpy(dtype=float)
    if v.size <= 1:
        return float("nan")
    return float(np.nanstd(v, ddof=1) / np.sqrt(v.size))


def _compute_wt_deltas(df: pd.DataFrame, metrics: list[str]) -> pd.DataFrame:
    out = df.copy()
    wt = out[out["mutation"] == "WT"].set_index(["structure", "replicate"])
    for metric in metrics:
        if metric not in out.columns or metric not in wt.columns:
            continue
        lookup = wt[metric]
        out[f"{metric}_delta"] = out.apply(
            lambda r: pd.to_numeric(r[metric], errors="coerce")
            - pd.to_numeric(lookup.get((r["structure"], r["replicate"]), np.nan), errors="coerce"),
            axis=1,
        )
    return out


def _aggregate_by_mutation(df: pd.DataFrame) -> pd.DataFrame:
    # Compute WT-deltas for structural ensemble metrics.
    df = _compute_wt_deltas(df, metrics=["contact_count", "hbond_count", "pocket_volume_proxy"])

    mut = df[df["mutation"] != "WT"].copy()
    if mut.empty:
        raise ValueError("No mutant rows found in ddg_full.csv (expected mutation != WT).")

    # Fold reduction should be constant per mutation; take first non-null.
    mut["fold_reduction"] = pd.to_numeric(mut["fold_reduction"], errors="coerce")
    if mut["fold_reduction"].isna().all():
        raise ValueError("fold_reduction missing for all mutants in ddg_full.csv.")

    agg = (
        mut.groupby("mutation", as_index=False)
        .agg(
            fold_reduction=("fold_reduction", "first"),
            n_replicates=("replicate", "nunique"),
            ddg=("ddg", "mean"),
            ddg_sem=("ddg", _sem),
            ddg_vdw=("ddg_vdw", "mean"),
            ddg_vdw_sem=("ddg_vdw", _sem),
            ddg_electrostatic=("ddg_electrostatic", "mean"),
            ddg_electrostatic_sem=("ddg_electrostatic", _sem),
            ddg_gb=("ddg_gb", "mean"),
            ddg_gb_sem=("ddg_gb", _sem),
            pocket_volume_proxy_delta=("pocket_volume_proxy_delta", "mean"),
            pocket_volume_proxy_delta_sem=("pocket_volume_proxy_delta", _sem),
            contact_count_delta=("contact_count_delta", "mean"),
            contact_count_delta_sem=("contact_count_delta", _sem),
            hbond_count_delta=("hbond_count_delta", "mean"),
            hbond_count_delta_sem=("hbond_count_delta", _sem),
        )
        .reset_index(drop=True)
    )

    agg["log10_fold_reduction"] = np.log10(agg["fold_reduction"].astype(float))
    agg["is_combo"] = agg["mutation"].astype(str).str.contains(r"\+")
    return agg


def _scatter(ax, df: pd.DataFrame, y: str, yerr: str, title: str, ylabel: str) -> None:
    # Separate singles vs combos for visual grouping.
    singles = df[~df["is_combo"]]
    combos = df[df["is_combo"]]

    def _plot(sub: pd.DataFrame, color: str, label: str, marker: str) -> None:
        ax.errorbar(
            sub["log10_fold_reduction"],
            sub[y],
            yerr=sub[yerr],
            fmt=marker,
            ms=6,
            lw=0,
            elinewidth=1,
            capsize=2,
            color=color,
            alpha=0.9,
            label=label,
        )

    _plot(singles, color="#1f77b4", label="Single", marker="o")
    _plot(combos, color="#d62728", label="Combination", marker="s")

    ax.set_title(title, fontsize=11, fontweight="bold")
    ax.set_xlabel(r"$\log_{10}(\mathrm{fold\ reduction})$")
    ax.set_ylabel(ylabel)
    ax.grid(True, alpha=0.25)


def plot_global_signatures(ddg_full_csv: Path, output_png: Path, annotate_top_n: int) -> None:
    import matplotlib.pyplot as plt

    df = pd.read_csv(ddg_full_csv)
    agg = _aggregate_by_mutation(df)

    # Order for annotation: highest fold reduction first.
    top = agg.sort_values("fold_reduction", ascending=False).head(max(0, int(annotate_top_n)))
    annotate_set = set(top["mutation"].astype(str).tolist())

    fig, axes = plt.subplots(2, 3, figsize=(14, 8), constrained_layout=True)
    specs = [
        ("ddg", "ddg_sem", "Total ΔΔG", "ΔΔG (MM/GBSA units)"),
        ("ddg_vdw", "ddg_vdw_sem", "van der Waals ΔΔG", "ΔΔG_vdw"),
        ("ddg_gb", "ddg_gb_sem", "Polar solvation ΔΔG", "ΔΔG_GB"),
        ("pocket_volume_proxy_delta", "pocket_volume_proxy_delta_sem", "Pocket volume shift", "Δ pocket volume (Å^3)"),
        ("contact_count_delta", "contact_count_delta_sem", "Contact network shift", "Δ contact count"),
        ("hbond_count_delta", "hbond_count_delta_sem", "H-bond network shift", "Δ H-bond count"),
    ]

    for ax, (y, yerr, title, ylabel) in zip(axes.ravel(), specs):
        _scatter(ax, agg, y=y, yerr=yerr, title=title, ylabel=ylabel)
        for _, row in agg.iterrows():
            mut = str(row["mutation"])
            if mut not in annotate_set:
                continue
            x = float(row["log10_fold_reduction"])
            yy = float(row[y])
            if not np.isfinite(x) or not np.isfinite(yy):
                continue
            ax.annotate(mut, (x, yy), xytext=(4, 4), textcoords="offset points", fontsize=8)

    # Legend in the first axis.
    axes[0, 0].legend(loc="best", fontsize=9, frameon=False)

    output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_png, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {output_png}")


def main() -> int:
    parser = argparse.ArgumentParser(description="Plot fold-reduction vs MM/GBSA/structural signatures.")
    parser.add_argument("--ddg-full", type=Path, default=Path("results/analysis/binding_energy/tables/ddg_full.csv"))
    parser.add_argument("--output", type=Path, default=Path("results/plots/manuscript_global_signatures.png"))
    parser.add_argument("--annotate-top-n", type=int, default=5, help="Annotate top-N highest fold-reduction mutations.")
    args = parser.parse_args()

    if not args.ddg_full.exists():
        raise FileNotFoundError(f"Missing: {args.ddg_full}")

    plot_global_signatures(args.ddg_full, args.output, annotate_top_n=args.annotate_top_n)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
