from __future__ import annotations

import pandas as pd


def plot_ddg_vs_fold_reduction(ddg_df: pd.DataFrame, paths) -> None:
    """Plot ΔΔG vs fold reduction for FEP results.

    Args:
        ddg_df: DataFrame with columns: mutation, ddg, fold_reduction, replicate.
        paths: Project paths object with plots directory.
    """
    import matplotlib.pyplot as plt
    import numpy as np

    mut_df = ddg_df[ddg_df["mutation"] != "WT"].dropna(subset=["ddg", "fold_reduction"])
    if mut_df.empty:
        return

    by_mutation = mut_df.groupby("mutation", as_index=False).agg(
        ddg_mean=("ddg", "mean"),
        ddg_sem=("ddg", lambda x: x.std() / np.sqrt(len(x)) if len(x) > 1 else 0),
        fold_reduction=("fold_reduction", "first"),
    )

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    ax = axes[0]
    ax.errorbar(
        by_mutation["fold_reduction"],
        by_mutation["ddg_mean"],
        yerr=by_mutation["ddg_sem"],
        fmt="o",
        color="#2a6f97",
        capsize=3,
        markersize=8,
        alpha=0.8,
    )
    for _, row in by_mutation.iterrows():
        ax.annotate(
            row["mutation"],
            (row["fold_reduction"], row["ddg_mean"]),
            textcoords="offset points",
            xytext=(5, 5),
            fontsize=7,
            alpha=0.7,
        )
    ax.set_xlabel("Fold Reduction (FC)")
    ax.set_ylabel("ΔΔG (kJ/mol)")
    ax.set_title("ΔΔG vs Fold Reduction")
    ax.axhline(0, color="#888888", linewidth=0.5, linestyle="--")

    if len(by_mutation) > 2:
        x = by_mutation["fold_reduction"].values
        y = by_mutation["ddg_mean"].values
        mask = np.isfinite(x) & np.isfinite(y)
        if mask.sum() > 2:
            coeffs = np.polyfit(x[mask], y[mask], 1)
            x_line = np.linspace(x[mask].min(), x[mask].max(), 100)
            y_line = np.polyval(coeffs, x_line)
            ax.plot(x_line, y_line, "--", color="#888888", alpha=0.6, label="Linear fit")

    ax = axes[1]
    log_fold = np.log10(by_mutation["fold_reduction"] + 1)
    ax.errorbar(
        log_fold,
        by_mutation["ddg_mean"],
        yerr=by_mutation["ddg_sem"],
        fmt="o",
        color="#7CB342",
        capsize=3,
        markersize=8,
        alpha=0.8,
    )
    for i, row in by_mutation.iterrows():
        ax.annotate(
            row["mutation"],
            (np.log10(row["fold_reduction"] + 1), row["ddg_mean"]),
            textcoords="offset points",
            xytext=(5, 5),
            fontsize=7,
            alpha=0.7,
        )
    ax.set_xlabel("log10(Fold Reduction + 1)")
    ax.set_ylabel("ΔΔG (kJ/mol)")
    ax.set_title("ΔΔG vs log10(Fold Reduction)")
    ax.axhline(0, color="#888888", linewidth=0.5, linestyle="--")

    if len(by_mutation) > 2:
        x = log_fold.values
        y = by_mutation["ddg_mean"].values
        mask = np.isfinite(x) & np.isfinite(y)
        if mask.sum() > 2:
            coeffs = np.polyfit(x[mask], y[mask], 1)
            x_line = np.linspace(x[mask].min(), x[mask].max(), 100)
            y_line = np.polyval(coeffs, x_line)
            ax.plot(x_line, y_line, "--", color="#888888", alpha=0.6, label="Linear fit")

    corr = by_mutation["ddg_mean"].corr(by_mutation["fold_reduction"])
    log_corr = by_mutation["ddg_mean"].corr(log_fold)
    fig.suptitle(
        f"DOR: ΔΔG Correlation with Resistance\n"
        f"(Pearson r={corr:.3f}, log-Pearson r={log_corr:.3f})",
        y=1.02,
    )

    fig.tight_layout()
    fig.savefig(paths.plots / "ddg_vs_fold_reduction.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_lambda_profiles(lambda_summary_df: pd.DataFrame, paths) -> None:
    """Plot cumulative ΔG profiles across lambda windows for sanity checks."""
    import matplotlib.pyplot as plt

    if lambda_summary_df.empty:
        return

    for leg in sorted(lambda_summary_df["leg"].dropna().unique()):
        leg_df = lambda_summary_df[lambda_summary_df["leg"] == leg].copy()
        if leg_df.empty:
            continue
        fig, ax = plt.subplots(figsize=(8, 5))
        for mutation in sorted(leg_df["mutation"].dropna().unique()):
            mdf = leg_df[leg_df["mutation"] == mutation].sort_values("window_index")
            ax.plot(
                mdf["window_index"],
                mdf["cumulative_delta_g_mean"],
                linewidth=1.5 if mutation == "WT" else 1.0,
                alpha=0.9 if mutation == "WT" else 0.45,
                label=mutation if mutation == "WT" else None,
                color="#1f77b4" if mutation == "WT" else "#999999",
            )
        ax.set_title(f"Cumulative ΔG Profile Across Lambda Windows ({leg} leg)")
        ax.set_xlabel("Window index")
        ax.set_ylabel("Cumulative ΔG (kJ/mol)")
        ax.axhline(0.0, color="#888888", linewidth=0.7, linestyle="--")
        if "WT" in set(leg_df["mutation"]):
            ax.legend(frameon=False)
        fig.tight_layout()
        fig.savefig(
            paths.plots / f"lambda_profile_{leg}.png",
            dpi=150,
            bbox_inches="tight",
        )
        plt.close(fig)
