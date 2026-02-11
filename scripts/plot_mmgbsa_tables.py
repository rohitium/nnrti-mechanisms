#!/usr/bin/env python3
from __future__ import annotations

import os
from pathlib import Path


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def main() -> int:
    root = _repo_root()
    os.environ.setdefault("MPLCONFIGDIR", str(root / ".mplconfig"))

    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt

    rep_path = root / "results" / "mmgbsa_replicate_metrics.csv"
    tbl_path = root / "results" / "table1_like_energy_components.csv"
    if not rep_path.exists():
        raise FileNotFoundError(rep_path)
    if not tbl_path.exists():
        raise FileNotFoundError(tbl_path)

    rep = pd.read_csv(rep_path)
    tbl = pd.read_csv(tbl_path)

    out_dir = root / "results" / "plots"
    out_dir.mkdir(parents=True, exist_ok=True)

    # --- Plot 1: replicate-level component bars (mean +- std over snapshots)
    comp_cols = [
        ("binding_dg_vdw", "vdW"),
        ("binding_dg_electrostatic", "Electrostatics"),
        ("binding_dg_gb", "GB (polar solvation)"),
        ("binding_dg_sa", "SA (nonpolar)"),
        ("binding_dg", "Total"),
    ]

    rep = rep.sort_values(["mutation", "replicate"]).reset_index(drop=True)
    mutations = rep["mutation"].tolist()
    snapshots = rep.get("mmgbsa_snapshots", pd.Series([np.nan] * len(rep))).dropna().unique()
    snapshots_note = ""
    if len(snapshots) == 1:
        snapshots_note = f" (n={int(snapshots[0])} snapshots)"
    elif len(snapshots) > 1:
        snapshots_note = f" (n snapshots: {', '.join(str(int(x)) for x in sorted(snapshots))})"

    x = np.arange(len(mutations))
    width = 0.16

    fig, ax = plt.subplots(figsize=(11, 4.8))
    for i, (col, label) in enumerate(comp_cols):
        y = rep[col].astype(float).to_numpy()
        yerr = rep.get(f"{col}_std", pd.Series([np.nan] * len(rep))).astype(float).to_numpy()
        ax.bar(x + (i - 2) * width, y, width=width, label=label, alpha=0.9)
        # Avoid matplotlib warnings for all-NaN yerr.
        if np.isfinite(yerr).any():
            ax.errorbar(x + (i - 2) * width, y, yerr=yerr, fmt="none", ecolor="black", elinewidth=1, capsize=2)

    ax.axhline(0.0, color="black", linewidth=1)
    ax.set_xticks(x, labels=[f"{m}" for m in mutations])
    ax.set_ylabel("Energy component (kJ/mol)")
    ax.set_title(f"MM/GBSA components{snapshots_note}")
    ax.legend(ncol=3, frameon=False)
    ax.grid(axis="y", linestyle=":", alpha=0.4)
    plt.tight_layout()
    fig.savefig(out_dir / "mmgbsa_components_by_mutation.png", dpi=300)
    plt.close(fig)

    # --- Plot 2: WT-referenced ddG bars (mutant - WT) with error bars.
    #
    # Prefer across-replicate SD when there are multiple replicates. If only one
    # replicate exists, fall back to propagated snapshot SEM:
    #   sem(ddG) = sqrt(sem(mut)^2 + sem(wt)^2)
    ddg_comps = [
        ("binding_dg_vdw", "vdW"),
        ("binding_dg_electrostatic", "Electrostatics"),
        ("binding_dg_gb", "GB (polar solvation)"),
        ("binding_dg_sa", "SA (nonpolar)"),
        ("binding_dg", "Total"),
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
            for col, _label in ddg_comps:
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
            **{f"ddg_{col}_mean": (f"ddg_{col}", "mean") for col, _ in ddg_comps},
            **{f"ddg_{col}_std": (f"ddg_{col}", "std") for col, _ in ddg_comps},
            **{f"ddg_{col}_sem_prop": (f"ddg_{col}_sem_prop", "mean") for col, _ in ddg_comps},
        )
        ddg_sum = ddg_sum.sort_values(["fold_reduction", "mutation"], ascending=[False, True]).reset_index(drop=True)

        mutations2 = ddg_sum["mutation"].tolist()
        x2 = np.arange(len(mutations2))
        fig, ax = plt.subplots(figsize=(11, 4.2))
        for i, (col, label) in enumerate(ddg_comps):
            y = ddg_sum[f"ddg_{col}_mean"].astype(float).to_numpy()
            # Error bar rule: SD across replicates if n>1 else propagated snapshot SEM.
            n_reps = ddg_sum["n_reps"].astype(int).to_numpy()
            yerr_sd = ddg_sum[f"ddg_{col}_std"].astype(float).to_numpy()
            yerr_sem = ddg_sum[f"ddg_{col}_sem_prop"].astype(float).to_numpy()
            yerr = np.where((n_reps > 1) & np.isfinite(yerr_sd), yerr_sd, yerr_sem)

            ax.bar(x2 + (i - 2) * width, y, width=width, label=label, alpha=0.9)
            if np.isfinite(yerr).any():
                ax.errorbar(
                    x2 + (i - 2) * width,
                    y,
                    yerr=yerr,
                    fmt="none",
                    ecolor="black",
                    elinewidth=1,
                    capsize=2,
                )
        ax.axhline(0.0, color="black", linewidth=1)
        ax.set_xticks(x2, labels=mutations2)
        ax.set_ylabel("ΔΔG component (kJ/mol) vs WT")
        ax.set_title(
            "∆∆G components (Y188L - WT)"
        )
        ax.legend(ncol=3, frameon=False)
        ax.grid(axis="y", linestyle=":", alpha=0.4)
        plt.tight_layout()
        fig.savefig(out_dir / "mmgbsa_ddg_components_vs_wt.png", dpi=300)
        plt.close(fig)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
