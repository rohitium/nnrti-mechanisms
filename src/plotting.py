from __future__ import annotations

import pandas as pd


def _safe_linear_fit(x_values, y_values):
    """Return (x_line, y_line, r) for a stable linear fit, or None."""
    import numpy as np

    x = np.asarray(x_values, dtype=float)
    y = np.asarray(y_values, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    if mask.sum() < 3:
        return None
    x = x[mask]
    y = y[mask]

    # Guard against degenerate/near-constant inputs that cause LAPACK/SVD failures.
    if np.nanstd(x) < 1e-12 or np.nanstd(y) < 1e-12:
        return None

    try:
        coeffs = np.polyfit(x, y, 1)
    except Exception:
        return None
    x_line = np.linspace(x.min(), x.max(), 100)
    y_line = np.polyval(coeffs, x_line)
    r = pd.Series(x).corr(pd.Series(y))
    return x_line, y_line, r


def cleanup_legacy_plots(paths) -> None:
    """Remove stale plot files from deprecated diagnostics."""
    legacy = [
        "ddg_vs_fold_reduction.png",
        "lambda_profile_complex.png",
        "lambda_profile_solvent.png",
        "lambda_hist_overlay_complex.png",
        "lambda_hist_overlay_solvent.png",
    ]
    for name in legacy:
        p = paths.plots / name
        if p.exists():
            p.unlink()


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

    fig, ax = plt.subplots(figsize=(13, 6))

    # Requested orientation: free energy on x-axis, susceptibility on y-axis.
    ax.errorbar(
        by_mutation["ddg_mean"],
        by_mutation["fold_reduction"],
        xerr=by_mutation["ddg_sem"],
        fmt="o",
        color="#2a6f97",
        capsize=3,
        markersize=8,
        alpha=0.85,
    )

    fit = _safe_linear_fit(
        by_mutation["ddg_mean"].values,
        by_mutation["fold_reduction"].values,
    )
    if fit is not None:
        x_line, y_line, _ = fit
        ax.plot(x_line, y_line, "--", color="#777777", alpha=0.7, label="Linear fit")

    # Spread labels with varied offsets to reduce overlap without extra deps.
    label_offsets = [
        (8, 8),
        (8, -10),
        (-10, 8),
        (-12, -10),
        (10, 16),
        (-10, 16),
        (12, -18),
        (-14, -18),
    ]
    sorted_idx = by_mutation.sort_values(["fold_reduction", "ddg_mean"]).index.tolist()
    for i, idx in enumerate(sorted_idx):
        row = by_mutation.loc[idx]
        dx, dy = label_offsets[i % len(label_offsets)]
        ax.annotate(
            row["mutation"],
            (row["ddg_mean"], row["fold_reduction"]),
            textcoords="offset points",
            xytext=(dx, dy),
            fontsize=8,
            alpha=0.9,
            bbox={
                "boxstyle": "round,pad=0.2",
                "facecolor": "white",
                "alpha": 0.7,
                "edgecolor": "none",
            },
        )

    ax.set_xlabel("ΔΔG (kJ/mol)")
    ax.set_ylabel("Fold Reduction (FC)")
    corr = by_mutation["ddg_mean"].corr(by_mutation["fold_reduction"])
    if pd.notna(corr):
        ax.set_title(
            f"DOR Resistance vs Binding Free-Energy Shift (Pearson r={corr:.3f})"
        )
    else:
        ax.set_title("DOR Resistance vs Binding Free-Energy Shift")
    ax.margins(x=0.12, y=0.12)
    if fit is not None:
        ax.legend(frameon=False, loc="upper left")

    fig.tight_layout()
    fig.savefig(paths.plots / "ddg_vs_fold_reduction.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


def _add_metric_delta_columns(ddg_df: pd.DataFrame) -> pd.DataFrame:
    df = ddg_df.copy()
    wt_df = df[df["mutation"] == "WT"].set_index(["structure", "replicate"])
    for metric in ("contact_count", "hbond_count", "pocket_volume_proxy"):
        if metric not in df.columns or metric not in wt_df.columns:
            continue
        wt_lookup = wt_df[metric]
        df[f"{metric}_delta"] = df.apply(
            lambda row: row[metric]
            - wt_lookup.get((row["structure"], row["replicate"]), float("nan")),
            axis=1,
        )
    return df


def plot_all_metrics_vs_fold_reduction(ddg_df: pd.DataFrame, paths) -> None:
    """Generate a multi-panel summary plot for all available metrics."""
    import matplotlib.pyplot as plt
    import numpy as np

    if ddg_df.empty:
        return

    df = _add_metric_delta_columns(ddg_df)
    mut_df = df[df["mutation"] != "WT"].copy()
    if mut_df.empty:
        return

    metric_specs = [
        ("ddg", "ΔΔG (kJ/mol)"),
        ("contact_count_delta", "Δ Contacts (mut - WT)"),
        ("hbond_count_delta", "Δ H-bonds (mut - WT)"),
        ("pocket_volume_proxy_delta", "Δ Pocket Volume Proxy (mut - WT)"),
    ]
    available = [m for m in metric_specs if m[0] in mut_df.columns]
    if not available:
        return

    agg = {"fold_reduction": ("fold_reduction", "first")}
    for col, _ in available:
        agg[col] = (col, "mean")
        agg[f"{col}_std"] = (col, "std")
        agg[f"{col}_n"] = (col, "count")
    by_mut = mut_df.groupby("mutation", as_index=False).agg(**agg)
    for col, _ in available:
        std_col = f"{col}_std"
        n_col = f"{col}_n"
        sem_col = f"{col}_sem"
        by_mut[sem_col] = by_mut[std_col] / np.sqrt(by_mut[n_col].clip(lower=1))
        by_mut[sem_col] = by_mut[sem_col].fillna(0.0)

    n = len(available)
    ncols = 2 if n > 1 else 1
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(7 * ncols, 4.8 * nrows))
    axes_arr = np.atleast_1d(axes).flatten()

    for i, (col, label) in enumerate(available):
        ax = axes_arr[i]
        sub = by_mut.dropna(subset=[col, "fold_reduction"])
        if sub.empty:
            ax.set_visible(False)
            continue
        sem_col = f"{col}_sem"
        ax.errorbar(
            sub["fold_reduction"],
            sub[col],
            yerr=sub[sem_col] if sem_col in sub.columns else None,
            fmt="o",
            color="#2a6f97",
            markersize=6,
            capsize=3,
            alpha=0.85,
        )
        for _, row in sub.iterrows():
            ax.annotate(
                row["mutation"],
                (row["fold_reduction"], row[col]),
                textcoords="offset points",
                xytext=(5, 4),
                fontsize=7,
                alpha=0.85,
            )
        fit = _safe_linear_fit(sub["fold_reduction"].values, sub[col].values)
        if fit is None:
            ax.set_title(f"{label} vs Fold Reduction")
        else:
            x_line, y_line, r = fit
            ax.plot(x_line, y_line, "--", color="#777777", alpha=0.65)
            if pd.notna(r):
                try:
                    from scipy import stats

                    _, pvalue = stats.pearsonr(
                        sub[col].values, sub["fold_reduction"].values
                    )
                except Exception:
                    pvalue = float("nan")
                r2 = float(r * r)
                if np.isfinite(pvalue):
                    ax.set_title(f"{label} vs Fold Reduction (R²={r2:.3f}, p={pvalue:.3g})")
                else:
                    ax.set_title(f"{label} vs Fold Reduction (R²={r2:.3f})")
            else:
                ax.set_title(f"{label} vs Fold Reduction")
        ax.set_xlabel("Fold Reduction (FC)")
        ax.set_ylabel(label)
        ax.margins(x=0.08, y=0.16)
        if col == "ddg":
            y = sub[col].to_numpy(dtype=float)
            y_min = float(np.nanmin(y))
            y_max = float(np.nanmax(y))
            y_span = y_max - y_min
            sem_max = (
                float(np.nanmax(sub[sem_col].to_numpy(dtype=float)))
                if sem_col in sub.columns
                else 0.0
            )
            pad = max(1.0, 0.6 * y_span, 3.0 * sem_max)
            ax.set_ylim(y_min - pad, y_max + pad)

    for j in range(len(available), len(axes_arr)):
        axes_arr[j].set_visible(False)

    fig.tight_layout()
    fig.savefig(
        paths.plots / "all_metrics_vs_fold_reduction.png",
        dpi=150,
        bbox_inches="tight",
    )
    plt.close(fig)


def plot_free_energy_convergence(convergence_df: pd.DataFrame, paths) -> None:
    """Plot WT convergence of cumulative ΔG vs fraction of collected samples."""
    import matplotlib.pyplot as plt

    if convergence_df.empty:
        return
    wt = convergence_df[convergence_df["mutation"] == "WT"].copy()
    if wt.empty:
        return

    summary = (
        wt.groupby(["leg", "sample_fraction"], as_index=False)
        .agg(
            delta_g_mean=("delta_g_kj_mol", "mean"),
            delta_g_std=("delta_g_kj_mol", "std"),
            n_replicates=("replicate", "nunique"),
        )
        .sort_values(["leg", "sample_fraction"])
    )
    if summary.empty:
        return

    fig, ax = plt.subplots(figsize=(8.5, 5))
    colors = {"complex": "#1f77b4", "solvent": "#ff7f0e"}
    for leg in sorted(summary["leg"].dropna().unique()):
        sub = summary[summary["leg"] == leg]
        ax.plot(
            sub["sample_fraction"],
            sub["delta_g_mean"],
            marker="o",
            linewidth=2.0,
            markersize=4,
            color=colors.get(leg, "#555555"),
            label=f"WT {leg}",
        )
        if "delta_g_std" in sub.columns:
            lo = sub["delta_g_mean"] - sub["delta_g_std"].fillna(0.0)
            hi = sub["delta_g_mean"] + sub["delta_g_std"].fillna(0.0)
            ax.fill_between(
                sub["sample_fraction"],
                lo,
                hi,
                color=colors.get(leg, "#555555"),
                alpha=0.14,
                linewidth=0,
            )

    ax.set_title("WT Free-Energy Convergence Across Sample Fraction")
    ax.set_xlabel("Fraction of production samples used")
    ax.set_ylabel("Estimated ΔG (kJ/mol)")
    ax.set_xlim(0.1, 1.0)
    ax.legend(frameon=False, loc="best")
    fig.tight_layout()
    fig.savefig(
        paths.plots / "free_energy_convergence_wt.png",
        dpi=150,
        bbox_inches="tight",
    )
    plt.close(fig)


def plot_boundness_qc(boundness_df: pd.DataFrame, paths) -> None:
    """Plot per-replicate minimum ligand-protein distance QC."""
    import matplotlib.pyplot as plt
    import numpy as np

    if boundness_df.empty:
        return
    df = boundness_df.copy().sort_values(["mutation", "replicate"])
    labels = [f"{m} r{int(r)}" for m, r in zip(df["mutation"], df["replicate"])]
    x = np.arange(len(df), dtype=float)
    threshold = float(df["bound_threshold_angstrom"].dropna().iloc[0]) if "bound_threshold_angstrom" in df.columns and df["bound_threshold_angstrom"].notna().any() else 6.0

    fig, ax = plt.subplots(figsize=(max(10, 0.45 * len(df)), 4.8))
    y_start = df["min_distance_start_angstrom"].to_numpy(dtype=float)
    ax.scatter(x, y_start, color="#d62728", s=34, label="Start structure min distance")

    if "min_distance_trajectory_angstrom" in df.columns:
        y_traj = df["min_distance_trajectory_angstrom"].to_numpy(dtype=float)
        mask = np.isfinite(y_traj)
        if np.any(mask):
            ax.scatter(
                x[mask],
                y_traj[mask],
                color="#1f77b4",
                s=34,
                marker="^",
                label="Trajectory sampled min distance",
            )

    for i, flag in enumerate(df.get("qc_flag", [])):
        if str(flag) != "OK":
            ax.annotate(
                str(flag),
                (x[i], y_start[i] if np.isfinite(y_start[i]) else threshold),
                textcoords="offset points",
                xytext=(0, 8),
                ha="center",
                fontsize=7,
                color="#aa0000",
            )

    ax.axhline(threshold, color="#444444", linestyle="--", linewidth=1.0)
    ax.text(
        0.99,
        0.98,
        f"Bound threshold = {threshold:.1f} Å",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=8,
        bbox={"boxstyle": "round,pad=0.2", "facecolor": "white", "alpha": 0.75, "edgecolor": "#666666"},
    )
    ax.set_title("Boundness QC: Minimum Ligand-Protein Distance per Replicate")
    ax.set_ylabel("Min distance (Å)")
    ax.set_xlabel("Mutation / replicate")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=70, ha="right", fontsize=7)
    ax.legend(frameon=False, loc="upper left")
    fig.tight_layout()
    fig.savefig(paths.plots / "boundness_qc_min_distance.png", dpi=160, bbox_inches="tight")
    plt.close(fig)


def plot_si_figure_s1_like(position_df: pd.DataFrame, paths) -> None:
    """S1-like mutation landscape: position frequency + resistance intensity."""
    import matplotlib.pyplot as plt

    if position_df.empty:
        return
    df = position_df.sort_values("position").copy()

    fig, axes = plt.subplots(2, 1, figsize=(10, 6), sharex=True)
    ax1, ax2 = axes
    ax1.bar(
        df["position"],
        df["n_mutations"],
        width=4.0,
        color="#4c78a8",
        alpha=0.85,
    )
    ax1.set_ylabel("Mutation count")
    ax1.set_title("Fig S1-like: Mutation Position Landscape (Our Dataset)")

    ax2.plot(
        df["position"],
        df["mean_log10_fold_plus1"],
        marker="o",
        linewidth=2.0,
        color="#f58518",
    )
    ax2.set_ylabel("Mean log10(Fold+1)")
    ax2.set_xlabel("RT position")
    ax2.grid(alpha=0.2, linestyle=":")

    fig.tight_layout()
    fig.savefig(paths.plots / "fig_s1_like_mutation_landscape.png", dpi=160, bbox_inches="tight")
    plt.close(fig)


def plot_si_figure_s2_like_rmsd(rmsd_df: pd.DataFrame, paths) -> None:
    """S2-like RMSD figure: C-alpha backbone RMSD vs simulation time."""
    import matplotlib.pyplot as plt
    import numpy as np

    if rmsd_df.empty:
        return
    df = rmsd_df.copy()
    if "time_ps" in df.columns and df["time_ps"].notna().any():
        x = df["time_ps"] / 1000.0
        xlab = "Time (ns)"
    else:
        x = df["frame_index"].astype(float)
        xlab = "Frame index"
    df["x"] = x

    muts = sorted(df["mutation"].dropna().unique().tolist())
    max_panels = min(4, len(muts))
    muts = muts[:max_panels]
    ncols = 2 if len(muts) > 1 else 1
    nrows = int(np.ceil(len(muts) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(11, 4.2 * nrows), sharex=False, sharey=True)
    axes_arr = np.atleast_1d(axes).flatten()

    for i, mut in enumerate(muts):
        ax = axes_arr[i]
        sub = df[df["mutation"] == mut].copy()
        for rep in sorted(sub["replicate"].dropna().unique()):
            rep_df = sub[sub["replicate"] == rep].sort_values("x")
            ax.plot(rep_df["x"], rep_df["ca_rmsd_angstrom"], linewidth=1.2, alpha=0.9, label=f"rep {int(rep)}")
        ax.set_title(str(mut))
        ax.set_ylabel("Cα RMSD (Å)")
        ax.set_xlabel(xlab)
        ax.grid(alpha=0.2, linestyle=":")
        ax.legend(frameon=False, fontsize=7)

    for j in range(len(muts), len(axes_arr)):
        axes_arr[j].set_visible(False)

    fig.suptitle("Fig S2-like: Backbone Cα RMSD Profiles", y=1.01)
    fig.tight_layout()
    fig.savefig(paths.plots / "fig_s2_like_ca_rmsd.png", dpi=160, bbox_inches="tight")
    plt.close(fig)


def plot_com_distance_convergence(com_df: pd.DataFrame, paths) -> None:
    """Plot DOR-RT COM distance convergence profiles by mutation/replicate."""
    import matplotlib.pyplot as plt
    import numpy as np

    if com_df.empty:
        return
    df = com_df.copy()
    if "time_ps" in df.columns and df["time_ps"].notna().any():
        df["x"] = df["time_ps"] / 1000.0
        xlab = "Time (ns)"
    else:
        df["x"] = df["frame_index"].astype(float)
        xlab = "Frame index"

    muts = sorted(df["mutation"].dropna().unique().tolist())
    max_panels = min(4, len(muts))
    muts = muts[:max_panels]
    ncols = 2 if len(muts) > 1 else 1
    nrows = int(np.ceil(len(muts) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(11, 4.2 * nrows), sharex=False, sharey=True)
    axes_arr = np.atleast_1d(axes).flatten()

    for i, mut in enumerate(muts):
        ax = axes_arr[i]
        sub = df[df["mutation"] == mut].copy()
        for rep in sorted(sub["replicate"].dropna().unique()):
            rep_df = sub[sub["replicate"] == rep].sort_values("x")
            ax.plot(rep_df["x"], rep_df["com_distance_angstrom"], linewidth=1.3, alpha=0.9, label=f"rep {int(rep)}")
        ax.set_title(str(mut))
        ax.set_ylabel("DOR-RT COM distance (Å)")
        ax.set_xlabel(xlab)
        ax.grid(alpha=0.2, linestyle=":")
        ax.legend(frameon=False, fontsize=7)

    for j in range(len(muts), len(axes_arr)):
        axes_arr[j].set_visible(False)

    fig.suptitle("DOR-RT Center-of-Mass Distance Convergence", y=1.01)
    fig.tight_layout()
    fig.savefig(paths.plots / "com_distance_convergence.png", dpi=160, bbox_inches="tight")
    plt.close(fig)


def plot_lambda_profiles(
    lambda_summary_df: pd.DataFrame,
    paths,
    significance_df: pd.DataFrame | None = None,
) -> None:
    """Plot WT cumulative λ profiles with replicate-derived SEM error bars."""
    import matplotlib.pyplot as plt
    import numpy as np

    if lambda_summary_df.empty:
        return

    wt = lambda_summary_df[lambda_summary_df["mutation"] == "WT"].copy()
    if wt.empty:
        return

    n_windows = int(wt["window_index"].max()) + 1
    if n_windows <= 0:
        return

    fig, axes = plt.subplots(2, 1, figsize=(9, 8), sharex=True)
    ax_total, ax_comp = axes
    leg_colors = {"complex": "#1f77b4", "solvent": "#2ca02c"}

    switch_x = min(4.0 / float(n_windows), 1.0)
    for ax in axes:
        ax.axvline(switch_x, color="#444444", linestyle=":", linewidth=1.0, alpha=0.9)
        ax.axvspan(0.0, switch_x, color="#4c78a8", alpha=0.06)
        ax.axvspan(switch_x, 1.0, color="#f58518", alpha=0.06)
        ax.axhline(0.0, color="#999999", linewidth=0.6, linestyle="--")

    for leg in sorted(wt["leg"].dropna().unique()):
        leg_df = wt[wt["leg"] == leg].sort_values("window_index").copy()
        x_progress = (leg_df["window_index"].to_numpy() + 1) / float(n_windows)
        x_plot = np.concatenate([[0.0], x_progress])

        nrep = leg_df["n_replicates"].to_numpy(dtype=float)
        nrep = np.where(nrep > 0, nrep, 1.0)

        total_mean = leg_df["cumulative_delta_g_mean"].to_numpy(dtype=float)
        total_sem = (
            leg_df["cumulative_delta_g_std"].fillna(0.0).to_numpy(dtype=float)
            / np.sqrt(nrep)
        )
        total_mean = np.concatenate([[0.0], total_mean])
        total_sem = np.concatenate([[0.0], total_sem])
        ax_total.errorbar(
            x_plot,
            total_mean,
            yerr=total_sem,
            color=leg_colors.get(leg, "#555555"),
            marker="o",
            markersize=3.5,
            linewidth=2.0,
            capsize=2,
            label=f"{leg} total",
        )

        e_mean = leg_df["cumulative_electrostatic_delta_g_mean"].to_numpy(dtype=float)
        e_sem = (
            leg_df["cumulative_electrostatic_delta_g_std"].fillna(0.0).to_numpy(dtype=float)
            / np.sqrt(nrep)
        )
        s_mean = leg_df["cumulative_steric_delta_g_mean"].to_numpy(dtype=float)
        s_sem = (
            leg_df["cumulative_steric_delta_g_std"].fillna(0.0).to_numpy(dtype=float)
            / np.sqrt(nrep)
        )
        e_mean = np.concatenate([[0.0], e_mean])
        e_sem = np.concatenate([[0.0], e_sem])
        s_mean = np.concatenate([[0.0], s_mean])
        s_sem = np.concatenate([[0.0], s_sem])

        color = leg_colors.get(leg, "#555555")
        ax_comp.plot(x_plot, e_mean, color=color, linewidth=1.8, linestyle="--", label=f"{leg} electrostatic")
        ax_comp.fill_between(x_plot, e_mean - e_sem, e_mean + e_sem, color=color, alpha=0.12, linewidth=0)
        ax_comp.plot(x_plot, s_mean, color=color, linewidth=1.8, linestyle="-.", label=f"{leg} steric")
        ax_comp.fill_between(x_plot, s_mean - s_sem, s_mean + s_sem, color=color, alpha=0.08, linewidth=0)

    ax_total.set_title("WT Cumulative ΔG Profiles (Complex vs RT-only/solvent)")
    ax_total.set_ylabel("Total Cumulative ΔG (kJ/mol)")
    ax_total.legend(frameon=False, loc="best")

    ax_comp.set_ylabel("Component Cumulative ΔG (kJ/mol)")
    ax_comp.set_xlabel("Alchemical progress λ (0 = DOR+RT, 1 = RT-only)")
    ax_comp.legend(frameon=False, loc="best", ncol=2, fontsize=8)
    ax_comp.set_xlim(0.0, 1.0)
    ax_comp.set_xticks([0.0, 0.25, 0.5, 0.75, 1.0])
    ax_comp.set_xticklabels(["0", "0.25", "0.5", "0.75", "1"])

    if significance_df is not None and not significance_df.empty:
        sig = significance_df.set_index("metric")
        parts = []
        for metric in ("total", "electrostatic", "steric"):
            if metric in sig.index:
                p = sig.loc[metric, "welch_pvalue"]
                if pd.notna(p):
                    parts.append(f"{metric} p={p:.3g}")
        if parts:
            ax_total.text(
                0.01,
                0.02,
                "WT complex vs solvent (Welch t-test): " + "; ".join(parts),
                transform=ax_total.transAxes,
                fontsize=8,
                ha="left",
                va="bottom",
                bbox={
                    "boxstyle": "round,pad=0.25",
                    "facecolor": "white",
                    "alpha": 0.78,
                    "edgecolor": "#666666",
                },
            )

    fig.tight_layout()
    fig.savefig(
        paths.plots / "lambda_profile_wt.png",
        dpi=150,
        bbox_inches="tight",
    )
    plt.close(fig)
