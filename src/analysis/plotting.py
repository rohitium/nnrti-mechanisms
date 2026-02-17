from __future__ import annotations

import pandas as pd

# 1 kcal = 4.184 kJ; OpenMM internally reports kJ/mol.
_KJ_TO_KCAL = 1.0 / 4.184


def cleanup_legacy_plots(paths) -> None:
    """Best-effort cleanup of deprecated plot artifacts."""
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

    def _load_structural_fallback() -> pd.DataFrame:
        candidates = [
            paths.results / "structural_metrics.csv",
            paths.results / ".checkpoints" / ".checkpoint_structural_metrics.csv",
            paths.results / ".checkpoint_structural_metrics.csv",
        ]
        for p in candidates:
            if not p.exists():
                continue
            try:
                sdf = pd.read_csv(p)
            except Exception:
                continue
            if sdf.empty:
                continue
            if "error" in sdf.columns:
                sdf = sdf[sdf["error"].isna()].copy()
            needed = {"structure", "mutation", "replicate", "fold_reduction"}
            if not needed.issubset(set(sdf.columns)):
                continue
            return _add_metric_delta_columns(sdf)
        return pd.DataFrame()

    if ddg_df.empty:
        df = _load_structural_fallback()
    else:
        df = ddg_df.copy()

    if df.empty:
        return

    mut_df = df[df["mutation"] != "WT"].copy()
    if mut_df.empty:
        return

    fold_valid = mut_df["fold_reduction"].notna().sum() if "fold_reduction" in mut_df.columns else 0
    if fold_valid < 2:
        fallback_df = _load_structural_fallback()
        if not fallback_df.empty:
            mut_df = fallback_df[fallback_df["mutation"] != "WT"].copy()

    energy_metric = "binding_dg" if "binding_dg" in mut_df.columns and mut_df["binding_dg"].notna().sum() >= 2 else "ddg"
    metric_specs = [
        (energy_metric, "Binding Energy (kcal/mol)"),
        ("pocket_volume_proxy", "Pocket Volume (A^3)"),
    ]

    available = []
    for col, label in metric_specs:
        if col not in mut_df.columns:
            continue
        valid = mut_df[[col, "fold_reduction"]].dropna()
        if len(valid) < 2:
            continue
        available.append((col, label))
    if not available:
        return

    agg = {"fold_reduction": ("fold_reduction", "first")}
    for col, _ in available:
        agg[col] = (col, "mean")
        agg[f"{col}_std"] = (col, "std")
        agg[f"{col}_n"] = (col, "count")

    by_mut = mut_df.groupby("mutation", as_index=False).agg(**agg)
    by_mut["fold_reduction"] = pd.to_numeric(by_mut["fold_reduction"], errors="coerce")
    by_mut["fold_reduction_log10"] = np.where(
        by_mut["fold_reduction"] > 0,
        np.log10(by_mut["fold_reduction"]),
        np.nan,
    )

    for energy_col in ("ddg", "binding_dg"):
        if energy_col in by_mut.columns:
            by_mut[energy_col] = by_mut[energy_col] * _KJ_TO_KCAL
            if f"{energy_col}_std" in by_mut.columns:
                by_mut[f"{energy_col}_std"] = by_mut[f"{energy_col}_std"] * _KJ_TO_KCAL

    wt_df = df[df["mutation"] == "WT"].copy()
    wt_summary: dict[str, tuple[float, float]] = {}
    for col, _ in available:
        if col not in wt_df.columns:
            continue
        wt_vals = pd.to_numeric(wt_df[col], errors="coerce").dropna().to_numpy(dtype=float)
        if wt_vals.size == 0:
            continue
        if col in {"ddg", "binding_dg"}:
            wt_vals = wt_vals * _KJ_TO_KCAL
        wt_mean = float(np.nanmean(wt_vals))
        wt_sem = float(np.nanstd(wt_vals, ddof=1) / np.sqrt(wt_vals.size)) if wt_vals.size > 1 else 0.0
        wt_summary[col] = (wt_mean, wt_sem)

    for col, _ in available:
        std_col = f"{col}_std"
        n_col = f"{col}_n"
        sem_col = f"{col}_sem"
        by_mut[sem_col] = by_mut[std_col] / np.sqrt(by_mut[n_col].clip(lower=1))
        by_mut[sem_col] = by_mut[sem_col].fillna(0.0)

    n = len(available)
    ncols = 2 if n > 1 else 1
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(8 * ncols, 5.5 * nrows))
    axes_arr = np.atleast_1d(axes).flatten()

    point_colors = {
        "single": "#1f77b4",
        "combo": "#d62728",
    }
    title_label_map = {
        "ddg": "Binding Energy",
        "binding_dg": "Binding Energy",
        "pocket_volume_proxy": "Pocket Volume",
    }

    def _annotate_non_overlapping(ax, sub_df: pd.DataFrame, x_col: str, y_col: str) -> None:
        fig = ax.figure
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        placed_boxes = []
        offsets = [
            (6, 6),
            (6, -8),
            (-8, 6),
            (-8, -8),
            (12, 0),
            (-12, 0),
            (0, 12),
            (0, -12),
            (14, 8),
            (-14, 8),
            (14, -8),
            (-14, -8),
        ]
        for _, row in sub_df.iterrows():
            for dx, dy in offsets:
                txt = ax.annotate(
                    row["mutation"],
                    (row[x_col], row[y_col]),
                    textcoords="offset points",
                    xytext=(dx, dy),
                    fontsize=8,
                    alpha=0.9,
                    bbox=dict(boxstyle="round,pad=0.25", facecolor="white", alpha=0.7, edgecolor="none"),
                )
                fig.canvas.draw()
                box = txt.get_window_extent(renderer=renderer).expanded(1.04, 1.12)
                if any(box.overlaps(prev) for prev in placed_boxes):
                    txt.remove()
                    continue
                placed_boxes.append(box)
                break

    for i, (col, label) in enumerate(available):
        ax = axes_arr[i]
        sub = by_mut.dropna(subset=[col, "fold_reduction_log10"])
        if sub.empty:
            ax.set_visible(False)
            continue

        sem_col = f"{col}_sem"
        singles = sub[~sub["mutation"].astype(str).str.contains(r"\+")].copy()
        combos = sub[sub["mutation"].astype(str).str.contains(r"\+")].copy()

        single_plot = None
        combo_plot = None
        if not singles.empty:
            single_err = singles[sem_col] if sem_col in singles.columns else None
            single_plot = ax.errorbar(
                singles["fold_reduction_log10"],
                singles[col],
                yerr=single_err,
                fmt="o",
                color=point_colors["single"],
                markersize=8,
                capsize=4,
                capthick=1.5,
                elinewidth=1.5,
                alpha=0.9,
                label="Single DRM",
            )
        if not combos.empty:
            combo_err = combos[sem_col] if sem_col in combos.columns else None
            combo_plot = ax.errorbar(
                combos["fold_reduction_log10"],
                combos[col],
                yerr=combo_err,
                fmt="s",
                color=point_colors["combo"],
                markersize=7,
                capsize=4,
                capthick=1.5,
                elinewidth=1.5,
                alpha=0.9,
                label="DRM Combination",
            )

        _annotate_non_overlapping(ax, sub, "fold_reduction_log10", col)

        wt_line = None
        if col in wt_summary:
            wt_mean, wt_sem = wt_summary[col]
            wt_line = ax.axhline(wt_mean, color="#000000", linestyle="-", linewidth=1.2, alpha=0.9, label="WT")
            if wt_sem > 0:
                ax.axhspan(wt_mean - wt_sem, wt_mean + wt_sem, color="#999999", alpha=0.15, linewidth=0)

        fit = _safe_linear_fit(sub["fold_reduction_log10"].values, sub[col].values)
        title_label = title_label_map.get(col, label)
        base_title = f"{title_label} vs " + r"$\log_{10}(\mathrm{Fold\ Change})$"
        trend_line = None
        trend_label = "Trend line"
        if fit is None:
            ax.set_title(base_title, fontsize=11, fontweight="normal")
        else:
            x_line, y_line, r = fit
            trend_line = ax.plot(x_line, y_line, "--", color="#777777", alpha=0.9, linewidth=1.8)[0]
            if pd.notna(r):
                try:
                    from scipy import stats

                    _, pvalue = stats.pearsonr(
                        sub[col].values, sub["fold_reduction_log10"].values
                    )
                except Exception:
                    pvalue = float("nan")
                r2 = float(r * r)
                if np.isfinite(pvalue):
                    trend_label = f"Trend line (R^2={r2:.3f}, p={pvalue:.3g})"
                else:
                    trend_label = f"Trend line (R^2={r2:.3f})"
        ax.set_title(base_title, fontsize=11, fontweight="normal")
        if trend_line is not None:
            trend_line.set_label(trend_label)

        ax.set_xlabel(r"$\log_{10}(\mathrm{Fold\ Change})$", fontsize=10)
        ax.set_ylabel(label, fontsize=10)
        ax.margins(x=0.08, y=0.16)
        ax.grid(alpha=0.3, linestyle=":", linewidth=0.8)

        legend_handles = []
        if single_plot is not None:
            legend_handles.append(single_plot)
        if combo_plot is not None:
            legend_handles.append(combo_plot)
        if wt_line is not None:
            legend_handles.append(wt_line)
        if trend_line is not None:
            legend_handles.append(trend_line)
        if legend_handles:
            ax.legend(handles=legend_handles, frameon=True, fontsize=8, loc="upper left")

        if i == 0 or col in {"ddg", "binding_dg"}:
            y = sub[col].to_numpy(dtype=float)
            yerr = (
                sub[sem_col].to_numpy(dtype=float)
                if sem_col in sub.columns
                else np.zeros_like(y)
            )
            y_low = float(np.nanmin(y - yerr))
            y_high = float(np.nanmax(y + yerr))
            y_span = y_high - y_low
            pad = max(0.05, 0.08 * y_span)
            ax.set_ylim(y_low - pad, y_high + pad)

    for j in range(len(available), len(axes_arr)):
        axes_arr[j].set_visible(False)

    fig.tight_layout()
    fig.savefig(
        paths.plots / "all_metrics_vs_fold_reduction.png",
        dpi=300,
        bbox_inches="tight",
    )
    plt.close(fig)


def plot_boundness_qc(boundness_df: pd.DataFrame, paths) -> None:
    """Plot mutation-level boundness QC summary with replicate error bars."""
    import matplotlib.pyplot as plt
    import numpy as np

    if boundness_df.empty:
        return
    df = boundness_df.copy()

    summary = (
        df.groupby("mutation", as_index=False)
        .agg(
            start_mean=("min_distance_start_angstrom", "mean"),
            start_std=("min_distance_start_angstrom", "std"),
            traj_mean=("min_distance_trajectory_angstrom", "mean"),
            traj_std=("min_distance_trajectory_angstrom", "std"),
            n_reps=("replicate", "nunique"),
        )
        .fillna({"start_std": 0.0, "traj_std": 0.0})
        .sort_values(
            "mutation",
            key=lambda s: s.map(
                lambda m: (0, str(m)) if str(m) == "WT" else (2, str(m)) if "+" in str(m) else (1, str(m))
            ),
        )
        .reset_index(drop=True)
    )
    if summary.empty:
        return

    n = summary["n_reps"].clip(lower=1).to_numpy(dtype=float)
    summary["start_sem"] = summary["start_std"].to_numpy(dtype=float) / np.sqrt(n)
    summary["traj_sem"] = summary["traj_std"].to_numpy(dtype=float) / np.sqrt(n)

    x = np.arange(len(summary), dtype=float)
    fig, ax = plt.subplots(figsize=(max(10, 0.7 * len(summary)), 5.2))
    ax.errorbar(
        x - 0.12,
        summary["start_mean"].to_numpy(dtype=float),
        yerr=summary["start_sem"].to_numpy(dtype=float),
        fmt="o",
        color="#d62728",
        markersize=5,
        capsize=3,
        label="Start structure min distance",
    )
    ax.errorbar(
        x + 0.12,
        summary["traj_mean"].to_numpy(dtype=float),
        yerr=summary["traj_sem"].to_numpy(dtype=float),
        fmt="^",
        color="#1f77b4",
        markersize=5,
        capsize=3,
        label="Trajectory sampled min distance",
    )

    ax.set_ylabel("Min distance (Å)")
    ax.set_xlabel("Mutation")
    ax.set_ylim(top=3.2)
    ax.set_xticks(x)
    ax.set_xticklabels(summary["mutation"].tolist(), rotation=45, ha="right", fontsize=8)
    ax.legend(frameon=False, loc="upper left")
    ax.grid(axis="y", alpha=0.25, linestyle=":")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(paths.plots / "boundness_qc_min_distance.png", dpi=160, bbox_inches="tight")
    plt.close(fig)


def plot_simulation_convergence(
    rmsd_df: pd.DataFrame,
    com_df: pd.DataFrame,
    paths,
) -> None:
    """Generate RMSD and COM distance convergence panels by mutation."""
    import matplotlib.pyplot as plt

    if rmsd_df.empty and com_df.empty:
        return

    all_muts: set[str] = set()
    if not rmsd_df.empty:
        all_muts |= set(rmsd_df["mutation"].dropna().unique())
    if not com_df.empty:
        all_muts |= set(com_df["mutation"].dropna().unique())
    if not all_muts:
        return

    def _mutation_sort_key(m: str) -> tuple[int, str]:
        if m == "WT":
            return (0, m)
        if "+" in m:
            return (2, m)
        return (1, m)

    muts = sorted(all_muts, key=_mutation_sort_key)

    max_cols = 4
    ncols = min(len(muts), max_cols)
    nrows = int((len(muts) + max_cols - 1) // max_cols)

    def _infer_total_ns_from_state_csv(safe_label: str, replicate: int) -> float | None:
        import pandas as pd

        state_csv = (
            paths.results
            / "md_runs"
            / str(safe_label)
            / f"rep_{int(replicate):02d}"
            / f"{safe_label}_rep{int(replicate):02d}_md_state.csv"
        )
        if not state_csv.exists():
            return None
        try:
            sdf = pd.read_csv(state_csv)
        except Exception:
            return None
        step_col = None
        for c in ('#"Step"', "Step"):
            if c in sdf.columns:
                step_col = c
                break
        if step_col is None or sdf.empty:
            return None
        max_step = pd.to_numeric(sdf[step_col], errors="coerce").dropna()
        if max_step.empty:
            return None
        return float(max_step.max()) * 2.0 / 1_000_000.0

    def _prep_x_fixed(df):
        import numpy as np
        import pandas as pd

        df = df.copy()
        if not {"safe_label", "replicate", "frame_index"}.issubset(df.columns):
            df["x"] = pd.to_numeric(df.get("time_ps", np.nan), errors="coerce") / 1000.0
            return df, "Time (ns)"

        parts = []
        for (safe_label, rep), sub in df.groupby(["safe_label", "replicate"], dropna=False):
            sub = sub.copy()
            max_frame = float(pd.to_numeric(sub["frame_index"], errors="coerce").max())
            total_ns = _infer_total_ns_from_state_csv(str(safe_label), int(rep))
            if total_ns is not None and np.isfinite(total_ns) and total_ns > 0 and max_frame > 0:
                sub["x"] = (pd.to_numeric(sub["frame_index"], errors="coerce") / max_frame) * total_ns
            else:
                t_ns = pd.to_numeric(sub.get("time_ps", np.nan), errors="coerce") / 1000.0
                # Legacy DCD metadata bug can inflate times by ~1000x.
                finite_max = float(np.nanmax(t_ns)) if np.isfinite(np.nanmax(t_ns)) else np.nan
                if np.isfinite(finite_max) and finite_max > 200.0:
                    t_ns = t_ns / 1000.0
                sub["x"] = t_ns
            parts.append(sub)
        df = pd.concat(parts, ignore_index=True) if parts else df
        return df, "Time (ns)"

    def _interp_mean_trace(sub: pd.DataFrame, y_col: str, n_grid: int = 200):
        import numpy as np

        reps = []
        for rep in sorted(sub["replicate"].dropna().unique()):
            rep_df = sub[sub["replicate"] == rep].sort_values("x")
            x = rep_df["x"].to_numpy(dtype=float)
            y = rep_df[y_col].to_numpy(dtype=float)
            if len(x) < 2:
                continue
            keep = np.r_[True, np.diff(x) > 0]
            x = x[keep]
            y = y[keep]
            if len(x) < 2:
                continue
            reps.append((x, y))

        if not reps:
            return None, None

        xmin = min(float(x.min()) for x, _ in reps)
        xmax = max(float(x.max()) for x, _ in reps)
        grid = np.linspace(xmin, xmax, n_grid)
        ys = []
        for x, y in reps:
            yi = np.interp(grid, x, y)
            yi[(grid < x.min()) | (grid > x.max())] = np.nan
            ys.append(yi)
        mean_y = np.nanmean(np.vstack(ys), axis=0)
        return grid, mean_y

    if not rmsd_df.empty:
        rdf, xlab = _prep_x_fixed(rmsd_df)
        fig_rmsd, axes_rmsd = plt.subplots(
            nrows=nrows,
            ncols=ncols,
            figsize=(5.5 * ncols, 3.5 * nrows),
            squeeze=False,
        )

        for plot_i, mut in enumerate(muts):
            row_i = plot_i // max_cols
            col_i = plot_i % max_cols
            ax = axes_rmsd[row_i, col_i]

            sub = rdf[rdf["mutation"] == mut]
            if sub.empty:
                ax.text(0.5, 0.5, "No data", ha="center", va="center", transform=ax.transAxes, color="#999999")
            else:
                color = "#1f77b4"
                for rep in sorted(sub["replicate"].dropna().unique()):
                    rep_df = sub[sub["replicate"] == rep].sort_values("x")
                    ax.plot(
                        rep_df["x"],
                        rep_df["ca_rmsd_angstrom"],
                        linewidth=1.0,
                        alpha=0.5,
                        color=color,
                    )
                x_mean, y_mean = _interp_mean_trace(sub, "ca_rmsd_angstrom")
                if x_mean is not None:
                    ax.plot(
                        x_mean,
                        y_mean,
                        linewidth=2.0,
                        alpha=1.0,
                        color=color,
                        label="replicate mean",
                    )

            ax.set_title(str(mut), fontsize=10, fontweight="bold")
            if col_i == 0:
                ax.set_ylabel("Cα RMSD (Å)", fontsize=9)
            ax.set_xlabel(xlab, fontsize=9)
            ax.grid(alpha=0.2, linestyle=":")

        for row_i in range(nrows):
            for col_i in range(ncols):
                idx = row_i * max_cols + col_i
                if idx >= len(muts):
                    axes_rmsd[row_i, col_i].set_visible(False)

        fig_rmsd.suptitle(
            "Root Mean-Squared Distance (RMSD) from the crystal structure",
            y=0.995,
            fontsize=14,
            fontweight="bold",
        )
        fig_rmsd.tight_layout()
        fig_rmsd.savefig(paths.plots / "rmsd_convergence.png", dpi=200, bbox_inches="tight")
        plt.close(fig_rmsd)

    if not com_df.empty:
        cdf, xlab = _prep_x_fixed(com_df)
        from matplotlib.ticker import FormatStrFormatter

        fig_com, axes_com = plt.subplots(
            nrows=nrows,
            ncols=ncols,
            figsize=(5.5 * ncols, 3.5 * nrows),
            squeeze=False,
        )

        for plot_i, mut in enumerate(muts):
            row_i = plot_i // max_cols
            col_i = plot_i % max_cols
            ax = axes_com[row_i, col_i]

            sub = cdf[cdf["mutation"] == mut]
            if sub.empty:
                ax.text(0.5, 0.5, "No data", ha="center", va="center", transform=ax.transAxes, color="#999999")
            else:
                color = "#d62728"
                for rep in sorted(sub["replicate"].dropna().unique()):
                    rep_df = sub[sub["replicate"] == rep].sort_values("x")
                    ax.plot(
                        rep_df["x"],
                        rep_df["com_distance_angstrom"],
                        linewidth=1.0,
                        alpha=0.5,
                        color=color,
                    )
                x_mean, y_mean = _interp_mean_trace(sub, "com_distance_angstrom")
                if x_mean is not None:
                    ax.plot(
                        x_mean,
                        y_mean,
                        linewidth=2.0,
                        alpha=1.0,
                        color=color,
                        label="replicate mean",
                    )

            ax.set_title(str(mut), fontsize=10, fontweight="bold")
            if col_i == 0:
                ax.set_ylabel("DOR–RT COM Distance (Å)", fontsize=9)
            ax.set_xlabel(xlab, fontsize=9)
            ax.yaxis.set_major_formatter(FormatStrFormatter("%.1f"))
            ax.grid(alpha=0.2, linestyle=":")

        for row_i in range(nrows):
            for col_i in range(ncols):
                idx = row_i * max_cols + col_i
                if idx >= len(muts):
                    axes_com[row_i, col_i].set_visible(False)

        fig_com.suptitle("DOR-RT Center of Mass (COM) Distances", y=0.995, fontsize=14, fontweight="bold")
        fig_com.tight_layout()
        fig_com.savefig(paths.plots / "com_distance_convergence.png", dpi=200, bbox_inches="tight")
        plt.close(fig_com)


def plot_all_mutation_pocket_volume_timeseries(pocket_df: pd.DataFrame, paths) -> None:
    """Plot WT-vs-mutant pocket-volume traces for each mutation."""
    import matplotlib.pyplot as plt
    import numpy as np

    if pocket_df.empty:
        return
    required = {"mutation", "replicate", "pocket_volume_proxy_angstrom3"}
    if not required.issubset(set(pocket_df.columns)):
        return

    def _infer_total_ns_from_state_csv(safe_label: str, replicate: int) -> float | None:
        state_csv = (
            paths.results
            / "md_runs"
            / str(safe_label)
            / f"rep_{int(replicate):02d}"
            / f"{safe_label}_rep{int(replicate):02d}_md_state.csv"
        )
        if not state_csv.exists():
            return None
        try:
            sdf = pd.read_csv(state_csv)
        except Exception:
            return None
        step_col = None
        for c in ('#"Step"', "Step"):
            if c in sdf.columns:
                step_col = c
                break
        if step_col is None or sdf.empty:
            return None
        max_step = pd.to_numeric(sdf[step_col], errors="coerce").dropna()
        if max_step.empty:
            return None
        return float(max_step.max()) * 2.0 / 1_000_000.0

    def _build_time_axis(df_in: pd.DataFrame) -> pd.DataFrame:
        out = df_in.copy()
        out["time_ns"] = pd.to_numeric(out.get("time_ps", np.nan), errors="coerce") / 1000.0
        if not {"safe_label", "replicate", "frame_index"}.issubset(out.columns):
            return out

        parts = []
        for (safe_label, rep), sub in out.groupby(["safe_label", "replicate"], dropna=False):
            sub = sub.copy()
            frame = pd.to_numeric(sub["frame_index"], errors="coerce")
            max_frame = float(frame.max()) if len(frame) else np.nan
            total_ns = _infer_total_ns_from_state_csv(str(safe_label), int(rep))
            if np.isfinite(max_frame) and max_frame > 0 and total_ns is not None and np.isfinite(total_ns) and total_ns > 0:
                sub["time_ns"] = (frame / max_frame) * total_ns
            else:
                t_ns = pd.to_numeric(sub.get("time_ps", np.nan), errors="coerce") / 1000.0
                finite_max = float(np.nanmax(t_ns)) if np.isfinite(np.nanmax(t_ns)) else np.nan
                if np.isfinite(finite_max) and finite_max > 200.0:
                    t_ns = t_ns / 1000.0
                sub["time_ns"] = t_ns.where(t_ns.notna(), frame * 0.01)
            parts.append(sub)
        return pd.concat(parts, ignore_index=True) if parts else out

    df = _build_time_axis(pocket_df)

    df["pocket_volume_nm3"] = pd.to_numeric(df["pocket_volume_proxy_angstrom3"], errors="coerce") * 1e-3
    df = df.dropna(subset=["mutation", "replicate", "time_ns", "pocket_volume_nm3"]).copy()
    if df.empty:
        return

    wt = df[df["mutation"] == "WT"].copy()
    muts = sorted([m for m in df["mutation"].unique() if m != "WT"], key=lambda m: (2, m) if "+" in str(m) else (1, m))
    if wt.empty or not muts:
        return

    out_dir = paths.plots / "pocket_volume_timeseries"
    out_dir.mkdir(parents=True, exist_ok=True)

    def _interp_mean_trace(sub: pd.DataFrame, n_grid: int = 220):
        reps = []
        for rep in sorted(sub["replicate"].dropna().unique()):
            rep_df = sub[sub["replicate"] == rep].sort_values("time_ns")
            x = rep_df["time_ns"].to_numpy(dtype=float)
            y = rep_df["pocket_volume_nm3"].to_numpy(dtype=float)
            if len(x) < 2:
                continue
            keep = np.r_[True, np.diff(x) > 0]
            x = x[keep]
            y = y[keep]
            if len(x) < 2:
                continue
            reps.append((x, y))
        if not reps:
            return None, None
        xmin = min(float(x.min()) for x, _ in reps)
        xmax = max(float(x.max()) for x, _ in reps)
        if not np.isfinite(xmin) or not np.isfinite(xmax) or xmax <= xmin:
            return None, None
        grid = np.linspace(xmin, xmax, n_grid)
        ys = []
        for x, y in reps:
            yi = np.interp(grid, x, y)
            yi[(grid < x.min()) | (grid > x.max())] = np.nan
            ys.append(yi)
        return grid, np.nanmean(np.vstack(ys), axis=0)

    for mut in muts:
        mdf = df[df["mutation"] == mut].copy()
        if mdf.empty:
            continue

        fig, ax = plt.subplots(figsize=(8.2, 5.0))

        for rep in sorted(wt["replicate"].dropna().unique()):
            rep_df = wt[wt["replicate"] == rep].sort_values("time_ns")
            ax.plot(rep_df["time_ns"], rep_df["pocket_volume_nm3"], color="#9a9a9a", alpha=0.35, linewidth=1.0)
        x_wt, y_wt = _interp_mean_trace(wt)
        if x_wt is not None:
            ax.plot(x_wt, y_wt, color="#4d4d4d", linewidth=2.2, label="WT mean")

        for rep in sorted(mdf["replicate"].dropna().unique()):
            rep_df = mdf[mdf["replicate"] == rep].sort_values("time_ns")
            ax.plot(rep_df["time_ns"], rep_df["pocket_volume_nm3"], color="#1f77b4", alpha=0.35, linewidth=1.0)
        x_mut, y_mut = _interp_mean_trace(mdf)
        if x_mut is not None:
            ax.plot(x_mut, y_mut, color="#1f77b4", linewidth=2.2, label=f"{mut} mean")

        ax.set_title(f"{mut} vs WT: Pocket Volume", fontsize=12, fontweight="bold")
        ax.set_xlabel("Time (ns)", fontsize=10)
        ax.set_ylabel("Pocket Volume (nm^3)", fontsize=10)
        ax.grid(alpha=0.25, linestyle=":")
        ax.legend(frameon=False, fontsize=9)
        fig.tight_layout()
        fig.savefig(out_dir / f"{str(mut).replace('+', '_')}_pocket_volume_timeseries.png", dpi=220, bbox_inches="tight")
        plt.close(fig)


def plot_si_figure_s1_like(pos_df: pd.DataFrame, paths) -> None:
    """Compatibility SI-style mutation landscape plot from position summary table."""
    import matplotlib.pyplot as plt
    import numpy as np

    if pos_df.empty:
        return

    # Accept multiple column naming conventions.
    pos_col = None
    for c in ("position", "residue_position", "site"):
        if c in pos_df.columns:
            pos_col = c
            break
    if pos_col is None:
        return

    fold_col = None
    for c in ("fold_reduction_mean", "mean_fold_reduction", "fold_reduction"):
        if c in pos_df.columns:
            fold_col = c
            break
    if fold_col is None:
        return

    n_col = None
    for c in ("n_mutations", "mutation_count", "n"):
        if c in pos_df.columns:
            n_col = c
            break

    df = pos_df.copy()
    df[pos_col] = pd.to_numeric(df[pos_col], errors="coerce")
    df[fold_col] = pd.to_numeric(df[fold_col], errors="coerce")
    df = df.dropna(subset=[pos_col, fold_col]).sort_values(pos_col)
    if df.empty:
        return

    x = df[pos_col].astype(int).to_numpy()
    y = np.log10(df[fold_col].clip(lower=1e-6).to_numpy())
    size = (
        40 + 40 * pd.to_numeric(df[n_col], errors="coerce").fillna(1).to_numpy()
        if n_col is not None
        else np.full(len(df), 80.0)
    )

    fig, ax = plt.subplots(figsize=(10, 4.8))
    ax.scatter(x, y, s=size, color="#1f77b4", alpha=0.85)
    ax.set_xlabel("RT residue position")
    ax.set_ylabel(r"$\log_{10}(\mathrm{FC})$")
    ax.set_title("Mutation Landscape by RT Position")
    ax.grid(alpha=0.25, linestyle=":")
    fig.tight_layout()
    fig.savefig(paths.plots / "fig_s1_like_mutation_landscape.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
