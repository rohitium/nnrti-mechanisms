from __future__ import annotations

import textwrap

import pandas as pd


def plot_delta_metrics(df: pd.DataFrame, paths) -> None:
    import matplotlib.pyplot as plt

    metric_labels = {
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
        fig, axes = plt.subplots(nrows=4, ncols=1, figsize=(fig_width, 8), sharex=True)
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
