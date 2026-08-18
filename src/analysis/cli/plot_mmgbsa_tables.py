#!/usr/bin/env python3
"""Plot MM/GBSA component and WT-referenced ddG panels.

Defaults follow the last-20-frame protocol (mmgbsa_snapshots=20), which is the
manuscript-facing analysis. The superseded 100-snapshot tables that used to sit at
results/mmgbsa_replicate_metrics.csv are archived under results/archive/.

Note: this script derives ddG in-script from the replicate table (mutant - WT per
replicate), so it takes no --ddg-csv; results/analysis/binding_energy/tables/ddg_full.csv
is not read here.
"""
from __future__ import annotations

import argparse
import os
from pathlib import Path

DEFAULT_REPLICATE_CSV = Path(
    "results/analysis/binding_energy/last20frames/mmgbsa_replicate_metrics_last20frames.csv"
)


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[3]


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot MM/GBSA component and WT-referenced ddG panels."
    )
    parser.add_argument(
        "--replicate-csv",
        type=Path,
        default=DEFAULT_REPLICATE_CSV,
        help="Replicate-level MM/GBSA metrics. Relative paths resolve against the repo root.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results/plots"),
        help="Directory for the generated PNGs. Relative paths resolve against the repo root.",
    )
    return parser.parse_args()


def main() -> int:
    args = _parse_args()
    root = _repo_root()
    os.environ.setdefault("MPLCONFIGDIR", str(root / ".mplconfig"))

    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt

    from src.analysis.units import KCAL_UNITS, read_energy_table

    rep_path = args.replicate_csv
    if not rep_path.is_absolute():
        rep_path = root / rep_path
    if not rep_path.exists():
        raise FileNotFoundError(rep_path)

    rep = read_energy_table(rep_path, KCAL_UNITS)

    out_dir = args.output_dir
    if not out_dir.is_absolute():
        out_dir = root / out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    # --- Plot 1: mutation-level component scatter with replicate SEM error bars.
    comp_cols = [
        ("binding_dg_vdw", "vdW", "#1f77b4"),
        ("binding_dg_electrostatic", "Electrostatics", "#2ca02c"),
        ("binding_dg_gb", "GB (polar solvation)", "#ff7f0e"),
        ("binding_dg_sa", "SA (nonpolar)", "#9467bd"),
        ("binding_dg", "Total", "#d62728"),
    ]

    # Order mutations: WT first, then single DRMs, then combinations.
    def _mutation_sort_key(m: str) -> tuple:
        if m == "WT":
            return (0, m)
        elif "+" in m:
            return (2, m)
        else:
            return (1, m)

    rep["_sort_key"] = rep["mutation"].apply(_mutation_sort_key)
    rep = rep.sort_values(["_sort_key", "mutation", "replicate"]).reset_index(drop=True)
    rep = rep.drop(columns=["_sort_key"])
    mut_order = sorted(rep["mutation"].unique(), key=_mutation_sort_key)

    # Aggregate replicate-level component means to mutation-level means ± SEM over replicates.
    agg = {"n_reps": ("replicate", "count")}
    for col, _, _ in comp_cols:
        agg[f"{col}_mean"] = (col, "mean")
        agg[f"{col}_std"] = (col, "std")
    by_mut = rep.groupby("mutation", as_index=False).agg(**agg)
    by_mut["mutation"] = pd.Categorical(by_mut["mutation"], categories=mut_order, ordered=True)
    by_mut = by_mut.sort_values("mutation").reset_index(drop=True)
    for col, _, _ in comp_cols:
        by_mut[f"{col}_sem"] = by_mut[f"{col}_std"] / np.sqrt(by_mut["n_reps"].clip(lower=1))
        by_mut[f"{col}_sem"] = by_mut[f"{col}_sem"].fillna(0.0)

    x = np.arange(len(by_mut), dtype=float)
    fig, axes = plt.subplots(
        nrows=len(comp_cols),
        ncols=1,
        figsize=(max(11, 0.75 * len(by_mut)), 2.4 * len(comp_cols) + 1.2),
        sharex=True,
        squeeze=False,
    )
    axes = axes[:, 0]
    for i, (col, label, color) in enumerate(comp_cols):
        ax = axes[i]
        y = by_mut[f"{col}_mean"].astype(float).to_numpy()
        yerr = by_mut[f"{col}_sem"].astype(float).to_numpy()
        ax.scatter(x, y, s=26, alpha=0.95, color=color, zorder=3)
        if np.isfinite(yerr).any():
            ax.errorbar(
                x, y, yerr=yerr, fmt="none", ecolor=color, elinewidth=1.1, capsize=2, alpha=0.95, zorder=2
            )
        ylo = y - np.where(np.isfinite(yerr), yerr, 0.0)
        yhi = y + np.where(np.isfinite(yerr), yerr, 0.0)
        finite = np.isfinite(ylo) & np.isfinite(yhi)
        if finite.any():
            dmin = float(np.nanmin(ylo[finite]))
            dmax = float(np.nanmax(yhi[finite]))
            span = max(1e-9, dmax - dmin)
            pad = 0.12 * span
            ax.set_ylim(dmin - pad, dmax + pad)
        wt_row = by_mut[by_mut["mutation"].astype(str) == "WT"]
        if not wt_row.empty and pd.notna(wt_row[f"{col}_mean"].iloc[0]):
            ax.axhline(
                float(wt_row[f"{col}_mean"].iloc[0]),
                color="#9a9a9a",
                linestyle="--",
                linewidth=1.0,
                alpha=0.8,
            )
        ax.set_ylabel(f"{label}\n(kcal/mol)", fontsize=9)
        ax.grid(axis="y", linestyle=":", alpha=0.35)

    axes[-1].set_xticks(x, labels=[str(m) for m in by_mut["mutation"]], rotation=45, ha="right", fontsize=8)
    axes[-1].set_xlabel("Mutation", fontsize=10)
    fig.suptitle("Binding Energy Components", fontsize=13, fontweight="bold", y=0.995)
    plt.tight_layout()
    fig.savefig(out_dir / "mmgbsa_components_by_mutation.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    # --- Plot 2: WT-referenced ddG scatter (mutant - WT) with SEM error bars.
    ddg_comps = [
        ("binding_dg_vdw", "vdW", "#1f77b4"),
        ("binding_dg_electrostatic", "Electrostatics", "#2ca02c"),
        ("binding_dg_gb", "GB (polar solvation)", "#ff7f0e"),
        ("binding_dg_sa", "SA (nonpolar)", "#9467bd"),
        ("binding_dg", "Total", "#d62728"),
    ]

    wt = rep[rep["mutation"] == "WT"].copy()
    muts = rep[rep["mutation"] != "WT"].copy()
    if not wt.empty and not muts.empty:
        wt_by_rep = {int(r): row for r, row in wt.set_index("replicate").iterrows()}
        wt_default = wt_by_rep.get(1)
        if wt_default is None:
            wt_default = next(iter(wt_by_rep.values()))

        rows = []
        for _, mrow in muts.iterrows():
            rep_id = int(mrow["replicate"])
            wrow = wt_by_rep.get(rep_id, wt_default)
            out = {"mutation": mrow["mutation"], "replicate": rep_id, "fold_reduction": mrow.get("fold_reduction", np.nan)}
            for col, _label, _color in ddg_comps:
                out[f"ddg_{col}"] = float(mrow[col]) - float(wrow[col])
                sem_m = float(mrow.get(f"{col}_sem", np.nan))
                sem_w = float(wrow.get(f"{col}_sem", np.nan))
                if np.isfinite(sem_m) and np.isfinite(sem_w):
                    out[f"ddg_{col}_sem_prop"] = float(np.sqrt(sem_m * sem_m + sem_w * sem_w))
                else:
                    out[f"ddg_{col}_sem_prop"] = np.nan
            rows.append(out)

        ddg = pd.DataFrame(rows)
        ddg_sum = ddg.groupby("mutation", as_index=False).agg(
            fold_reduction=("fold_reduction", "max"),
            n_reps=("replicate", "count"),
            **{f"ddg_{col}_mean": (f"ddg_{col}", "mean") for col, _, _ in ddg_comps},
            **{f"ddg_{col}_std": (f"ddg_{col}", "std") for col, _, _ in ddg_comps},
            **{f"ddg_{col}_sem_prop": (f"ddg_{col}_sem_prop", "mean") for col, _, _ in ddg_comps},
        )
        for col, _, _ in ddg_comps:
            ddg_sum[f"ddg_{col}_sem_rep"] = ddg_sum[f"ddg_{col}_std"] / np.sqrt(
                ddg_sum["n_reps"].clip(lower=1).astype(float)
            )
        ddg_sum = ddg_sum.sort_values(["fold_reduction", "mutation"], ascending=[False, True]).reset_index(drop=True)

        mutations2 = ddg_sum["mutation"].tolist()
        x2 = np.arange(len(mutations2))
        fig, axes = plt.subplots(
            nrows=len(ddg_comps),
            ncols=1,
            figsize=(max(11, 0.75 * len(mutations2)), 2.4 * len(ddg_comps) + 1.2),
            sharex=True,
            squeeze=False,
        )
        axes = axes[:, 0]
        for i, (col, label, color) in enumerate(ddg_comps):
            ax = axes[i]
            y = ddg_sum[f"ddg_{col}_mean"].astype(float).to_numpy()
            # Error bar rule: SEM across replicates if n>1 else propagated snapshot SEM.
            n_reps = ddg_sum["n_reps"].astype(int).to_numpy()
            yerr_sem_rep = ddg_sum[f"ddg_{col}_sem_rep"].astype(float).to_numpy()
            yerr_sem = ddg_sum[f"ddg_{col}_sem_prop"].astype(float).to_numpy()
            yerr = np.where((n_reps > 1) & np.isfinite(yerr_sem_rep), yerr_sem_rep, yerr_sem)
            ax.scatter(x2, y, s=26, alpha=0.95, color=color, zorder=3)
            if np.isfinite(yerr).any():
                ax.errorbar(
                    x2,
                    y,
                    yerr=yerr,
                    fmt="none",
                    ecolor=color,
                    elinewidth=1.1,
                    capsize=2,
                    alpha=0.95,
                    zorder=2,
                )
            ax.axhline(0.0, color="#9a9a9a", linestyle="--", linewidth=1.0, alpha=0.8)
            ax.set_ylabel(f"{label}\n(kcal/mol)", fontsize=9)
            ax.grid(axis="y", linestyle=":", alpha=0.35)
        axes[-1].set_xticks(x2, labels=mutations2, rotation=45, ha="right", fontsize=9)
        axes[-1].set_xlabel("Mutation", fontsize=10)
        fig.suptitle("Binding Energy Components (Mut - WT)", fontsize=13, fontweight="bold", y=0.995)
        plt.tight_layout()
        fig.savefig(out_dir / "mmgbsa_ddg_components_vs_wt.png", dpi=300, bbox_inches="tight")
        plt.close(fig)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
