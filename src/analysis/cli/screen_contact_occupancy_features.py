#!/usr/bin/env python3
"""Build residue-level occupancy_mean features and rank them against DOR susceptibility."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

from ..susceptibility import load_dor_susceptibilities
from .plot_dor_susceptibility_bars import (
    CATEGORY_COLORS,
    NEGATIVE_CONTROLS,
    POSITIVE_CONTROLS,
    UNCERTAIN_LIMITED,
)
from .plot_triplet_contact_story import _load_replicate_meta


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Rank residue-DOR occupancy_mean features by susceptibility correlation.")
    parser.add_argument("--manifest", type=Path, default=Path("manifests/md_manifest.csv"))
    parser.add_argument("--susceptibility-xlsx", type=Path, default=Path("data/DRM-susceptibilities.csv.xlsx"))
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results/analysis/contact_occupancy_feature_screen"),
    )
    parser.add_argument("--ligand-resname", type=str, default="2KW")
    parser.add_argument("--resid-offset", type=int, default=-3, help="auth = traj - resid_offset")
    parser.add_argument("--contact-cutoff", type=float, default=4.0)
    parser.add_argument("--window-ns", type=float, default=100.0)
    parser.add_argument(
        "--min-any-occupancy-mean",
        type=float,
        default=0.5,
        help="Keep residue features whose mutation-level occupancy_mean exceeds this threshold in at least one mutation.",
    )
    parser.add_argument(
        "--top-k",
        type=int,
        default=0,
        help="Number of top features to plot. Use 0 to plot all selected features.",
    )
    parser.add_argument(
        "--max-panels-per-figure",
        type=int,
        default=6,
        help="Maximum number of feature panels to place in a single figure before splitting into parts.",
    )
    parser.add_argument(
        "--reuse-existing-tables",
        action="store_true",
        help="Reuse cached occupancy feature tables in the output directory and only rerender the plots.",
    )
    return parser.parse_args()


def _fold_map_with_wt(xlsx_path: Path) -> dict[str, float]:
    df = load_dor_susceptibilities(xlsx_path)
    out = {str(row["mutation"]): float(row["dor_fold_reduction"]) for _, row in df.iterrows()}
    out["WT"] = 1.0
    return out


def _category_for_mutation(label: str) -> str:
    mutation = str(label).strip().upper()
    if mutation == "WT":
        return "WT"
    if mutation in NEGATIVE_CONTROLS:
        return "Negative control"
    if mutation in POSITIVE_CONTROLS:
        return "Positive control"
    if mutation in UNCERTAIN_LIMITED:
        return "Uncertain/limited data"
    return "Other"


def _category_color(label: str) -> str:
    cat = _category_for_mutation(label)
    if cat == "WT":
        return "#333333"
    if cat in CATEGORY_COLORS:
        return CATEGORY_COLORS[cat]
    return "#777777"


def _feature_label(resname: str, auth_resid: int) -> str:
    name = str(resname).strip().upper() or "UNK"
    return f"{name}{int(auth_resid)}"


def _wt_resname_map_from_manifest(manifest_csv: Path, resid_offset: int, auth_resids: Iterable[int] | None = None) -> dict[int, str]:
    repo_root = Path(__file__).resolve().parents[3]
    mf = pd.read_csv(manifest_csv)
    wt_rows = mf[mf["mutation"].astype(str).str.upper() == "WT"].sort_values("replicate", kind="stable")
    if wt_rows.empty:
        return {}
    output_json = wt_rows.iloc[0]["output_json"]
    out_json = Path(str(output_json))
    if not out_json.exists():
        marker = "nnrti-mechanisms/"
        text = str(out_json)
        if marker in text:
            out_json = repo_root / text.split(marker, 1)[1]
        else:
            out_json = repo_root / str(out_json)
    if not out_json.exists():
        return {}
    data = json.loads(out_json.read_text())
    topo = Path(str(data.get("analysis_topology_pdb", "")))
    if not topo.exists():
        marker = "nnrti-mechanisms/"
        text = str(topo)
        if marker in text:
            topo = repo_root / text.split(marker, 1)[1]
        else:
            topo = repo_root / str(topo)
    if not topo.exists():
        return {}

    wanted = None if auth_resids is None else {int(v) for v in auth_resids}
    mapping: dict[int, str] = {}
    for line in topo.read_text().splitlines():
        if not line.startswith("ATOM"):
            continue
        chain = line[21].strip()
        if chain != "A":
            continue
        resname = line[17:20].strip().upper()
        try:
            traj_resid = int(line[22:26].strip())
        except ValueError:
            continue
        auth_resid = int(traj_resid) - int(resid_offset)
        if wanted is not None and auth_resid not in wanted:
            continue
        mapping.setdefault(auth_resid, resname)
    return mapping


def _compute_contact_occupancy_tables(
    metas,
    ligand_resname: str,
    resid_offset: int,
    contact_cutoff: float,
    window_ns: float,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    import MDAnalysis as mda
    from MDAnalysis import transformations as trans
    from MDAnalysis.lib.distances import capped_distance

    rep_counts: dict[tuple[str, int, int], int] = {}
    rep_frames: dict[tuple[str, int], int] = {}
    res_meta: dict[int, str] = {}
    timing_rows: list[dict[str, object]] = []

    for m in sorted(metas, key=lambda x: (x.mutation, x.replicate)):
        u = mda.Universe(str(m.topology_pdb), str(m.analysis_dcd), format="DCD")
        lig = u.select_atoms(f"resname {ligand_resname} and not name H*")
        prot = u.select_atoms("protein and not name H*")
        if lig.n_atoms == 0 or prot.n_atoms == 0:
            continue

        try:
            anchor = u.select_atoms("protein")
            if anchor.n_atoms == 0:
                anchor = u.atoms
            u.trajectory.add_transformations(
                trans.NoJump(check_continuity=False),
                trans.center_in_box(anchor, center="geometry", wrap=False),
            )
        except Exception:
            pass

        n_frames = len(u.trajectory)
        if n_frames < 2:
            continue
        total_ns = float(m.total_ns) if np.isfinite(m.total_ns) and m.total_ns > 0 else float(window_ns)
        t_ns = np.linspace(0.0, total_ns, n_frames)
        keep_idx = np.where(t_ns <= float(window_ns))[0]
        if keep_idx.size < 2:
            keep_idx = np.arange(n_frames, dtype=int)

        atom_to_resid = np.asarray([int(a.resid) for a in prot.atoms], dtype=int)
        rep_key = (m.mutation, int(m.replicate))
        rep_frames[rep_key] = int(len(keep_idx))
        timing_rows.append(
            {
                "mutation": m.mutation,
                "replicate": int(m.replicate),
                "n_frames_total": int(n_frames),
                "n_frames_window": int(len(keep_idx)),
                "total_ns_used": float(total_ns),
                "timing_source": m.timing_source,
                "analysis_dcd": str(m.analysis_dcd),
            }
        )

        for fi in keep_idx.tolist():
            u.trajectory[int(fi)]
            pairs = capped_distance(
                prot.positions,
                lig.positions,
                max_cutoff=float(contact_cutoff),
                min_cutoff=0.0,
                box=u.dimensions,
                return_distances=False,
            )
            if pairs is None or len(pairs) == 0:
                continue
            touched = set(atom_to_resid[np.asarray(pairs)[:, 0]].tolist())
            for tr in touched:
                k = (m.mutation, int(m.replicate), int(tr))
                rep_counts[k] = rep_counts.get(k, 0) + 1
                if int(tr) not in res_meta:
                    ag = u.select_atoms(f"protein and resid {int(tr)}")
                    res_meta[int(tr)] = str(ag.residues[0].resname) if ag.n_atoms > 0 and len(ag.residues) else ""

    rep_rows: list[dict[str, object]] = []
    for (mutation, replicate), nfr in sorted(rep_frames.items()):
        touched_resids = sorted({k[2] for k in rep_counts.keys() if k[0] == mutation and k[1] == replicate})
        for tr in touched_resids:
            cnt = int(rep_counts.get((mutation, replicate, tr), 0))
            auth = int(tr) - int(resid_offset)
            rep_rows.append(
                {
                    "mutation": mutation,
                    "replicate": int(replicate),
                    "traj_resid": int(tr),
                    "auth_resid": int(auth),
                    "resname": str(res_meta.get(int(tr), "")),
                    "n_contact_frames": cnt,
                    "n_total_frames": int(nfr),
                    "occupancy": float(cnt / max(1, nfr)),
                }
            )
    rep_occ = pd.DataFrame(rep_rows)
    if rep_occ.empty:
        raise ValueError("No residue contacts were detected for the selected simulations.")

    mut_occ = (
        rep_occ.groupby(["mutation", "traj_resid", "auth_resid", "resname"], as_index=False)
        .agg(
            occupancy_mean=("occupancy", "mean"),
            occupancy_std=("occupancy", "std"),
            n_replicates=("replicate", "nunique"),
            n_contact_frames_total=("n_contact_frames", "sum"),
            n_total_frames_total=("n_total_frames", "sum"),
        )
    )
    mut_occ["occupancy_sem"] = mut_occ["occupancy_std"] / np.sqrt(mut_occ["n_replicates"].clip(lower=1))
    mut_occ["occupancy_pooled"] = (
        pd.to_numeric(mut_occ["n_contact_frames_total"], errors="coerce")
        / pd.to_numeric(mut_occ["n_total_frames_total"], errors="coerce").clip(lower=1)
    )
    timing_df = pd.DataFrame(timing_rows)
    return rep_occ, mut_occ, timing_df


def _rank_features(
    feature_df: pd.DataFrame,
    feature_cols: list[str],
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for col in feature_cols:
        sub = feature_df[["mutation", "dor_fold_reduction", "log10_dor_fold_reduction", col]].copy()
        sub[col] = pd.to_numeric(sub[col], errors="coerce")
        sub = sub[np.isfinite(sub["dor_fold_reduction"]) & np.isfinite(sub[col])]
        if len(sub) < 3:
            continue
        x_raw = sub["dor_fold_reduction"].to_numpy(dtype=float)
        x_log = sub["log10_dor_fold_reduction"].to_numpy(dtype=float)
        y = sub[col].to_numpy(dtype=float)
        pearson_raw = pearsonr(x_raw, y)
        pearson_log = pearsonr(x_log, y)
        spearman_raw = spearmanr(x_raw, y)
        rows.append(
            {
                "feature": col,
                "n": int(len(sub)),
                "pearson_r_fold": float(pearson_raw.statistic),
                "pearson_p_fold": float(pearson_raw.pvalue),
                "pearson_r_log10_fold": float(pearson_log.statistic),
                "pearson_p_log10_fold": float(pearson_log.pvalue),
                "spearman_rho_fold": float(spearman_raw.statistic),
                "spearman_p_fold": float(spearman_raw.pvalue),
                "abs_pearson_r_fold": float(abs(pearson_raw.statistic)),
                "feature_mean": float(np.nanmean(y)),
                "feature_std": float(np.nanstd(y, ddof=1)) if len(y) > 1 else np.nan,
            }
        )
    out = pd.DataFrame(rows)
    if out.empty:
        raise ValueError("No rankable occupancy_mean features were found.")
    return out.sort_values(["abs_pearson_r_fold", "pearson_r_fold"], ascending=[False, False], kind="stable").reset_index(drop=True)


def _plot_top_features(
    feature_df: pd.DataFrame,
    top_features: pd.DataFrame,
    sem_map: dict[str, pd.Series],
    output_png: Path,
    max_panels_per_figure: int,
) -> None:
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    def choose_offsets(df: pd.DataFrame) -> list[tuple[float, float]]:
        candidates = [
            (10.0, 8.0),
            (10.0, -8.0),
            (-10.0, 8.0),
            (-10.0, -8.0),
            (14.0, 0.0),
            (-14.0, 0.0),
            (6.0, 13.0),
            (6.0, -13.0),
        ]
        xlog = np.log10(df["dor_fold_reduction"].to_numpy(dtype=float))
        yvals = df["occupancy_mean"].to_numpy(dtype=float)
        xspan = max(float(np.nanmax(xlog) - np.nanmin(xlog)), 0.2)
        yspan = max(float(np.nanmax(yvals) - np.nanmin(yvals)), 0.1)
        anchors: list[tuple[float, float]] = []
        offsets: list[tuple[float, float]] = []
        for px, py in zip(xlog, yvals):
            best = candidates[0]
            best_score = -np.inf
            for dx, dy in candidates:
                ax = px + dx / 90.0
                ay = py + dy / 180.0
                if not anchors:
                    score = 0.0
                else:
                    d2 = [((ax - ox) / xspan) ** 2 + ((ay - oy) / max(yspan, 1e-6)) ** 2 for ox, oy in anchors]
                    score = min(d2)
                if score > best_score:
                    best_score = score
                    best = (dx, dy)
            anchors.append((px + best[0] / 90.0, py + best[1] / 180.0))
            offsets.append(best)
        return offsets

    top_features = top_features.reset_index(drop=True).copy()
    n_total = len(top_features)
    n_per_fig = max(1, int(max_panels_per_figure))
    n_figs = int(np.ceil(n_total / n_per_fig))
    output_png.parent.mkdir(parents=True, exist_ok=True)

    legend_handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor="#333333", markeredgecolor="#333333", label="WT"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor=CATEGORY_COLORS["Negative control"], markeredgecolor=CATEGORY_COLORS["Negative control"], label="Negative control"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor=CATEGORY_COLORS["Positive control"], markeredgecolor=CATEGORY_COLORS["Positive control"], label="Positive control"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor=CATEGORY_COLORS["Uncertain/limited data"], markeredgecolor=CATEGORY_COLORS["Uncertain/limited data"], label="Uncertain/limited data"),
    ]

    for fig_idx in range(n_figs):
        chunk = top_features.iloc[fig_idx * n_per_fig : (fig_idx + 1) * n_per_fig].copy()
        n = len(chunk)
        ncols = 2
        nrows = int(np.ceil(n / ncols))
        fig, axes = plt.subplots(nrows, ncols, figsize=(13.8, 3.9 * nrows), constrained_layout=True)
        axes = np.atleast_1d(axes).ravel()

        for ax, (_, row) in zip(axes, chunk.iterrows()):
            feature = str(row["feature"])
            df = feature_df[["mutation", "dor_fold_reduction", feature]].copy()
            df["occupancy_mean"] = pd.to_numeric(df[feature], errors="coerce")
            df["occupancy_sem"] = pd.to_numeric(sem_map[feature], errors="coerce")
            df = df[np.isfinite(df["dor_fold_reduction"]) & np.isfinite(df["occupancy_mean"])].copy()
            df = df.sort_values("dor_fold_reduction", kind="stable")
            offsets = choose_offsets(df)

            for (_, pt), (dx, dy) in zip(df.iterrows(), offsets):
                mutation = str(pt["mutation"])
                color = _category_color(mutation)
                px = float(pt["dor_fold_reduction"])
                py = float(pt["occupancy_mean"])
                pyerr = float(pt["occupancy_sem"]) if np.isfinite(pt["occupancy_sem"]) and float(pt["occupancy_sem"]) > 0.0 else None
                if pyerr is not None:
                    ax.errorbar(
                        px,
                        py,
                        yerr=pyerr,
                        fmt="none",
                        ecolor="#222222",
                        elinewidth=1.8,
                        capsize=4.0,
                        capthick=1.6,
                        alpha=0.9,
                        zorder=2,
                    )
                ax.scatter(
                    [px],
                    [py],
                    s=38,
                    facecolor=color,
                    edgecolor="white",
                    linewidth=0.6,
                    alpha=0.98,
                    zorder=3,
                )
                ax.annotate(
                    mutation,
                    xy=(px, py),
                    xytext=(dx, dy),
                    textcoords="offset points",
                    fontsize=7.0,
                    color=color,
                    va="center",
                    ha="left" if dx >= 0 else "right",
                    bbox={"boxstyle": "round,pad=0.12", "fc": "white", "ec": "none", "alpha": 0.72},
                    annotation_clip=True,
                )

            x = df["dor_fold_reduction"].to_numpy(dtype=float)
            y = df["occupancy_mean"].to_numpy(dtype=float)
            xlog = np.log10(x)
            if len(xlog) >= 2 and np.nanstd(xlog) > 0 and np.nanstd(y) > 0:
                coeffs = np.polyfit(xlog, y, 1)
                xline = np.logspace(float(np.nanmin(xlog)), float(np.nanmax(xlog)), 200)
                yline = coeffs[0] * np.log10(xline) + coeffs[1]
                ax.plot(xline, yline, color="#444444", linewidth=1.4, linestyle="--", alpha=0.9, zorder=2)

            ax.set_title(
                f"{feature.replace('occupancy_mean_', '')}\n"
                f"Pearson r(log10 fold) = {float(row['pearson_r_log10_fold']):.3f}, "
                f"R² = {float(row['pearson_r_log10_fold']) ** 2:.3f}, p = {float(row['pearson_p_log10_fold']):.3g}",
                fontsize=11,
            )
            ax.set_xlabel("Fold-change")
            ax.set_xscale("log")
            xticks = [1, 2, 5, 10, 20, 50, 100]
            x_present = [tick for tick in xticks if float(np.nanmin(x)) <= tick <= float(np.nanmax(x))]
            if 1.0 not in x_present and float(np.nanmin(x)) <= 1.0 <= float(np.nanmax(x)):
                x_present = [1.0] + x_present
            if x_present:
                ax.set_xticks(x_present)
                ax.set_xticklabels([f"{tick:g}" for tick in x_present])
            ax.set_ylabel("Occupancy Mean ± SEM")
            ax.set_xlim(left=0.8)
            ylo = float(np.nanmin(np.minimum(df["occupancy_mean"].to_numpy(dtype=float), (df["occupancy_mean"] - df["occupancy_sem"].fillna(0.0)).to_numpy(dtype=float))))
            yhi = float(np.nanmax(np.maximum(df["occupancy_mean"].to_numpy(dtype=float), (df["occupancy_mean"] + df["occupancy_sem"].fillna(0.0)).to_numpy(dtype=float))))
            span = max(yhi - ylo, 0.08)
            pad = 0.14 * span
            lower = max(-0.02, ylo - pad)
            upper = min(1.02, yhi + pad)
            if upper - lower < 0.08:
                mid = 0.5 * (upper + lower)
                lower = max(-0.02, mid - 0.04)
                upper = min(1.02, mid + 0.04)
            ax.set_ylim(lower, upper)
            ax.grid(alpha=0.22, linestyle=":")
            for spine in ("top", "right"):
                ax.spines[spine].set_visible(False)

        for ax in axes[n:]:
            ax.axis("off")

        axes[0].legend(handles=legend_handles, loc="lower right", frameon=True, fontsize=8)

        if n_figs == 1:
            png_path = output_png
            pdf_path = output_png.with_suffix(".pdf")
        else:
            stem = output_png.stem
            png_path = output_png.with_name(f"{stem}_part_{fig_idx + 1:02d}.png")
            pdf_path = output_png.with_name(f"{stem}_part_{fig_idx + 1:02d}.pdf")
        fig.savefig(png_path, dpi=300)
        fig.savefig(pdf_path)
        plt.close(fig)


def main() -> int:
    args = _parse_args()
    if not args.manifest.exists():
        raise FileNotFoundError(args.manifest)
    if not args.susceptibility_xlsx.exists():
        raise FileNotFoundError(args.susceptibility_xlsx)

    out_tables = args.output_dir / "tables"
    out_plots = args.output_dir / "plots"
    out_config = args.output_dir / "config"
    out_tables.mkdir(parents=True, exist_ok=True)
    out_plots.mkdir(parents=True, exist_ok=True)
    out_config.mkdir(parents=True, exist_ok=True)

    if bool(args.reuse_existing_tables):
        feature_matrix = pd.read_csv(out_tables / "occupancy_mean_feature_matrix.csv")
        sem_matrix = pd.read_csv(out_tables / "occupancy_mean_sem_matrix.csv")
        ranking = pd.read_csv(out_tables / "occupancy_feature_ranking.csv")
    else:
        fold_map = _fold_map_with_wt(args.susceptibility_xlsx)
        needed_mutations = set(fold_map.keys())
        metas = _load_replicate_meta(args.manifest, needed_mutations=needed_mutations)
        if not metas:
            raise ValueError("No replicate metadata found for fold-labeled mutations.")

        rep_occ, mut_occ, timing_df = _compute_contact_occupancy_tables(
            metas=metas,
            ligand_resname=str(args.ligand_resname),
            resid_offset=int(args.resid_offset),
            contact_cutoff=float(args.contact_cutoff),
            window_ns=float(args.window_ns),
        )

        rep_occ.to_csv(out_tables / "replicate_contact_occupancy.csv", index=False)
        mut_occ.to_csv(out_tables / "mutation_contact_occupancy.csv", index=False)
        timing_df.to_csv(out_tables / "timing_audit.csv", index=False)

        selected = (
            mut_occ.groupby(["traj_resid", "auth_resid", "resname"], as_index=False)["occupancy_mean"]
            .max()
            .rename(columns={"occupancy_mean": "max_occupancy_mean"})
        )
        selected = selected[selected["max_occupancy_mean"] > float(args.min_any_occupancy_mean)].copy()
        wt_name_map = _wt_resname_map_from_manifest(
            args.manifest,
            resid_offset=int(args.resid_offset),
            auth_resids=selected["auth_resid"].astype(int).tolist(),
        )
        selected["feature"] = selected.apply(
            lambda r: f"occupancy_mean_{_feature_label(str(wt_name_map.get(int(r['auth_resid']), r['resname'])), int(r['auth_resid']))}",
            axis=1,
        )
        selected = selected.sort_values(["auth_resid", "traj_resid"], kind="stable").reset_index(drop=True)
        selected.to_csv(out_tables / "selected_occupancy_features.csv", index=False)

        mut_occ = mut_occ.copy()
        mut_occ["feature"] = mut_occ.apply(
            lambda r: f"occupancy_mean_{_feature_label(str(wt_name_map.get(int(r['auth_resid']), r['resname'])), int(r['auth_resid']))}",
            axis=1,
        )
        use = mut_occ[mut_occ["feature"].isin(set(selected["feature"]))].copy()

        feature_matrix = (
            use.pivot_table(index="mutation", columns="feature", values="occupancy_mean", aggfunc="mean")
            .reset_index()
            .copy()
        )
        sem_matrix = (
            use.pivot_table(index="mutation", columns="feature", values="occupancy_sem", aggfunc="mean")
            .reset_index()
            .copy()
        )
        feature_matrix["dor_fold_reduction"] = feature_matrix["mutation"].map(fold_map)
        feature_matrix["log10_dor_fold_reduction"] = np.log10(feature_matrix["dor_fold_reduction"].astype(float))
        feature_matrix["category"] = feature_matrix["mutation"].map(_category_for_mutation)
        feature_matrix = feature_matrix.sort_values("dor_fold_reduction", kind="stable").reset_index(drop=True)
        feature_matrix.to_csv(out_tables / "occupancy_mean_feature_matrix.csv", index=False)
        sem_matrix.to_csv(out_tables / "occupancy_mean_sem_matrix.csv", index=False)

        feature_cols = [c for c in feature_matrix.columns if c.startswith("occupancy_mean_")]
        ranking = _rank_features(feature_matrix, feature_cols)
        ranking.to_csv(out_tables / "occupancy_feature_ranking.csv", index=False)

    feature_cols = [c for c in feature_matrix.columns if c.startswith("occupancy_mean_")]

    if int(args.top_k) <= 0 or int(args.top_k) >= len(ranking):
        top_features = ranking.copy()
        plot_stem = "all_occupancy_mean_features_vs_fold_change"
    else:
        top_features = ranking.head(int(args.top_k)).copy()
        plot_stem = f"top_{int(args.top_k)}_occupancy_mean_features_vs_fold_change"
    top_features.to_csv(out_tables / "top_occupancy_features.csv", index=False)

    sem_map = {
        col: sem_matrix.set_index("mutation")[col]
        for col in feature_cols
        if col in sem_matrix.columns
    }
    _plot_top_features(
        feature_df=feature_matrix,
        top_features=top_features,
        sem_map=sem_map,
        output_png=out_plots / f"{plot_stem}.png",
        max_panels_per_figure=int(args.max_panels_per_figure),
    )

    (out_config / "run_config.json").write_text(
        json.dumps(
            {
                "manifest": str(args.manifest),
                "susceptibility_xlsx": str(args.susceptibility_xlsx),
                "output_dir": str(args.output_dir),
                "ligand_resname": str(args.ligand_resname),
                "resid_offset": int(args.resid_offset),
                "contact_cutoff": float(args.contact_cutoff),
                "window_ns": float(args.window_ns),
                "min_any_occupancy_mean": float(args.min_any_occupancy_mean),
                "top_k": int(args.top_k),
                "max_panels_per_figure": int(args.max_panels_per_figure),
                "feature_definition": "mutation-level occupancy_mean = unweighted mean of per-replicate residue-DOR contact occupancies; contact defined as any heavy-atom pair within cutoff",
                "feature_selection_rule": "keep residues whose mutation-level occupancy_mean exceeds threshold in at least one mutation",
                "ranking_metric": "absolute Pearson correlation with raw DOR fold-change",
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
