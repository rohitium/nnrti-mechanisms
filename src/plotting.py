from __future__ import annotations

import textwrap

import pandas as pd


def plot_delta_metrics(df: pd.DataFrame, paths) -> None:
    import matplotlib.pyplot as plt

    metric_labels = {
        "binding_delta_g_kj_mol": "Δ Binding Free Energy (kJ/mol)",
        "binding_proxy_kj_mol": "Δ Binding Energy (kJ/mol)",
        "contact_count": "Δ Contacts",
        "hbond_count": "Δ H-bonds",
        "pocket_volume_proxy": "Δ Binding Pocket\nVolume (A^3)",
    }
    bar_colors = {"DOR": "#7CB342", "RPV": "#FFB74D"}

    def _wrap_xtick_label(label: str) -> str:
        label = label.replace("+", "+\n")
        return textwrap.fill(label, width=10)

    for structure in df["structure"].unique():
        struct_df = df[df["structure"] == structure].copy()
        struct_df["mutation_label"] = struct_df["mutation"].astype(str)
        mutation_order = (
            struct_df[struct_df["state"] == "MUT"][
                ["mutation_label", "mutation_order"]
            ]
            .drop_duplicates()
            .sort_values("mutation_order")
        )
        mutations = mutation_order["mutation_label"].tolist()
        fig_width = max(6, 0.45 * max(1, len(mutations)))
        fig, axes = plt.subplots(
            nrows=len(metric_labels), ncols=1, figsize=(fig_width, 10), sharex=True
        )
        bar_color = bar_colors.get(structure, "#2a6f97")
        for idx, (metric, label) in enumerate(metric_labels.items()):
            metric_df = struct_df[struct_df["metric"] == metric]
            if "replicate" in metric_df.columns:
                pivot = metric_df.pivot_table(
                    index=["mutation_label", "replicate"],
                    columns="state",
                    values="value",
                    aggfunc="first",
                )
                pivot = pivot.dropna(subset=["WT", "MUT"], how="any")
                pivot["delta"] = pivot["MUT"] - pivot["WT"]
                delta_stats = (
                    pivot["delta"]
                    .groupby(level=0)
                    .agg(["mean", "std", "count"])
                    .rename(columns={"mean": "delta_mean", "std": "delta_std"})
                )
                if mutations:
                    delta_stats = delta_stats.reindex(mutations)
                y = delta_stats["delta_mean"].values
                yerr = delta_stats["delta_std"].fillna(0.0).values
                axes[idx].bar(
                    range(len(mutations)),
                    y,
                    color=bar_color,
                    yerr=yerr,
                    capsize=2,
                )
            else:
                pivot = metric_df.pivot_table(
                    index="mutation_label",
                    columns="state",
                    values="value",
                    aggfunc="first",
                )
                if mutations:
                    pivot = pivot.reindex(mutations)
                delta = pivot["MUT"] - pivot["WT"]
                axes[idx].bar(range(len(mutations)), delta.values, color=bar_color)
            axes[idx].axhline(0, color="#444444", linewidth=0.8)
            axes[idx].set_ylabel(label)
            wrapped_labels = [_wrap_xtick_label(label) for label in mutations]
            axes[idx].set_xticks(range(len(mutations)))
            axes[idx].set_xticklabels(
                wrapped_labels, fontsize=6, rotation=45, ha="center"
            )
        axes[-1].set_xlabel("")
        fig.suptitle(f"{structure} (Δ = Mutation - WT)", y=0.98)
        fig.tight_layout(rect=[0, 0, 1, 0.96])
        fig.savefig(paths.plots / f"{structure.lower()}_delta_metrics.png", dpi=150)
        plt.close(fig)


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
