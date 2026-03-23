#!/usr/bin/env python3
"""Plot feature distributions for triplets matched to contact-story comparisons."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde

from ..susceptibility import load_dor_susceptibilities


def _parse_triplets(path: Path) -> list[tuple[str, str, str]]:
    triplets: list[tuple[str, str, str]] = []
    for line in path.read_text().splitlines():
        text = str(line).strip()
        if not text:
            continue
        toks = [tok.strip() for tok in text.replace(",", "|").split("|") if tok.strip()]
        if len(toks) != 3:
            raise ValueError(f"Triplet line must contain 3 mutations: {text}")
        triplets.append((toks[0], toks[1], toks[2]))
    if not triplets:
        raise ValueError(f"No triplets found in {path}")
    return triplets


def _fold_map(xlsx_path: Path) -> dict[str, float]:
    df = load_dor_susceptibilities(xlsx_path)
    out = {str(row["mutation"]): float(row["dor_fold_reduction"]) for _, row in df.iterrows()}
    out["WT"] = 1.0
    return out


def _label(mutation: str, fold_map: dict[str, float]) -> str:
    fold = float(fold_map.get(str(mutation), np.nan))
    if np.isfinite(fold):
        return f"{mutation} ({fold:.1f}x)"
    return str(mutation)


def _safe_name(mutation: str) -> str:
    return str(mutation).replace("+", "_").replace(" ", "_")


def _feature_slug(feature: str) -> str:
    text = str(feature).strip()
    if text.endswith("_angstrom"):
        text = text[: -len("_angstrom")]
    return text


def _feature_axis_label(feature: str) -> str:
    labels = {
        "ligand_palm_distance_angstrom": "Ligand Distance To NNIBP Palm Centroid (A)",
        "ligand_pose_rmsd_angstrom": "Ligand Pose RMSD To 4NCG (A)",
        "ligand_entrance_distance_angstrom": "Ligand Distance To NNIBP Entrance Centroid (A)",
        "ligand_pocket_center_distance_angstrom": "Ligand Distance To NNIBP Pocket Center (A)",
        "ligand_palm_depth_projection_angstrom": "Ligand Projection Along Entrance-To-Palm Axis (A)",
    }
    text = str(feature).strip()
    if text.startswith("residue_min_distance_") and text.endswith("_angstrom"):
        core = text[len("residue_min_distance_") : -len("_angstrom")]
        split = 0
        for idx, ch in enumerate(core):
            if ch.isdigit():
                split = idx
                break
        if split > 0:
            residue_name = core[:split]
            residue_id = core[split:]
            return f"Minimum Heavy-Atom Distance: {residue_name}{residue_id} To DOR (A)"
    return labels.get(str(feature), str(feature).replace("_", " "))


def _feature_title(feature: str) -> str:
    titles = {
        "ligand_palm_distance_angstrom": "Palm-Distance Histograms With KDE In Triplet Susceptibility Context",
        "ligand_pose_rmsd_angstrom": "Ligand Pose RMSD Histograms With KDE In Triplet Susceptibility Context",
        "ligand_entrance_distance_angstrom": "Entrance-Distance Histograms With KDE In Triplet Susceptibility Context",
        "ligand_pocket_center_distance_angstrom": "Pocket-Center Distance Histograms With KDE In Triplet Susceptibility Context",
        "ligand_palm_depth_projection_angstrom": "Palm-Depth Projection Histograms With KDE In Triplet Susceptibility Context",
    }
    text = str(feature).strip()
    if text.startswith("residue_min_distance_") and text.endswith("_angstrom"):
        core = text[len("residue_min_distance_") : -len("_angstrom")]
        split = 0
        for idx, ch in enumerate(core):
            if ch.isdigit():
                split = idx
                break
        if split > 0:
            residue_name = core[:split]
            residue_id = core[split:]
            return f"{residue_name}{residue_id}-To-DOR Distance Histograms With KDE In Triplet Susceptibility Context"
    return titles.get(str(feature), f"{str(feature).replace('_', ' ')} distributions")


def _mutation_summaries(frame_df: pd.DataFrame, feature: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    rep = (
        frame_df.groupby(["mutation", "replicate"], as_index=False)[feature]
        .agg(["count", "mean", "std", "median", "min", "max"])
        .reset_index()
    )
    rep = rep.rename(columns={"count": "n_frames", "std": f"{feature}_std", "mean": feature, "median": f"{feature}_median"})
    mut = (
        rep.groupby("mutation", as_index=False)[feature]
        .agg(["count", "mean", "std", "median", "min", "max"])
        .reset_index()
    )
    mut = mut.rename(columns={"count": "n_replicates", "std": f"{feature}_repstd", "mean": feature, "median": f"{feature}_median"})
    return rep, mut


def _plot_triplet(
    frame_df: pd.DataFrame,
    rep_df: pd.DataFrame,
    triplet: tuple[str, str, str],
    *,
    feature: str,
    fold_map: dict[str, float],
    output_png: Path,
) -> dict[str, object]:
    colors = {
        triplet[0]: "#1d3557",
        triplet[1]: "#6c757d",
        triplet[2]: "#d62828",
    }
    order = list(triplet)

    values = [frame_df.loc[frame_df["mutation"] == mut, feature].dropna().to_numpy(dtype=float) for mut in order]
    if any(v.size < 2 for v in values):
        raise ValueError(f"Insufficient feature samples for triplet {triplet}")

    lo = min(float(np.min(v)) for v in values)
    hi = max(float(np.max(v)) for v in values)
    pad = max(0.25, 0.08 * (hi - lo))
    grid = np.linspace(lo - pad, hi + pad, 300)
    bins = np.linspace(lo - pad, hi + pad, 24)

    fig, axes = plt.subplots(3, 1, figsize=(8.8, 6.6), sharex=True)
    if not isinstance(axes, np.ndarray):
        axes = np.asarray([axes])

    summary_rows: list[dict[str, object]] = []
    for ax, mut, pooled in zip(axes, order, values):
        rep = rep_df[rep_df["mutation"] == mut].copy().sort_values("replicate")
        ax.hist(
            pooled,
            bins=bins,
            density=True,
            color=colors[mut],
            alpha=0.28,
            edgecolor="white",
            linewidth=0.6,
        )
        kde = gaussian_kde(pooled)
        ax.plot(grid, kde(grid), color=colors[mut], linewidth=2.0)
        pooled_mean = float(np.mean(pooled))
        pooled_median = float(np.median(pooled))
        ax.axvline(pooled_mean, color=colors[mut], linestyle="--", linewidth=1.3, alpha=0.95)
        ax.axvline(pooled_median, color="#222222", linestyle=":", linewidth=1.1, alpha=0.9)
        ax.set_ylabel("Density")
        ax.set_title(_label(mut, fold_map), loc="left", fontsize=11)
        ax.grid(axis="y", alpha=0.22)
        summary_rows.append(
            {
                "triplet": "|".join(triplet),
                "mutation": mut,
                "n_frames": int(len(pooled)),
                "n_replicates": int(rep["replicate"].nunique()),
                "pooled_mean": pooled_mean,
                "pooled_median": pooled_median,
                "pooled_std": float(np.std(pooled, ddof=1)) if len(pooled) > 1 else 0.0,
                "pooled_q25": float(np.quantile(pooled, 0.25)),
                "pooled_q75": float(np.quantile(pooled, 0.75)),
                "replicate_mean_std": float(pd.to_numeric(rep[feature], errors="coerce").std(ddof=1)) if len(rep) > 1 else 0.0,
                "dor_fold_reduction": float(fold_map.get(mut, np.nan)),
            }
        )

    axes[-1].set_xlabel(_feature_axis_label(feature))
    fig.suptitle(_feature_title(feature), y=0.98)
    output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.97))
    fig.savefig(output_png, dpi=220, bbox_inches="tight")
    plt.close(fig)
    return {"plot_path": str(output_png), "summary_rows": summary_rows}


def main() -> int:
    parser = argparse.ArgumentParser(description="Plot feature distributions for contact-story triplets.")
    parser.add_argument(
        "--frame-feature-csv",
        type=Path,
        default=Path("results/analysis/ligand_pocket_features/tables/frame_features.csv"),
    )
    parser.add_argument(
        "--triplets-txt",
        type=Path,
        default=Path("results/analysis/triplet_contact_story_100ns/config/triplets.txt"),
    )
    parser.add_argument(
        "--susceptibility-xlsx",
        type=Path,
        default=Path("data/DRM-susceptibilities.csv.xlsx"),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results/analysis/palm_distance_triplets"),
    )
    parser.add_argument(
        "--feature-column",
        type=str,
        default="ligand_palm_distance_angstrom",
    )
    args = parser.parse_args()

    if not args.frame_feature_csv.exists():
        raise FileNotFoundError(args.frame_feature_csv)
    if not args.triplets_txt.exists():
        raise FileNotFoundError(args.triplets_txt)
    if not args.susceptibility_xlsx.exists():
        raise FileNotFoundError(args.susceptibility_xlsx)

    frame_df = pd.read_csv(args.frame_feature_csv)
    if args.feature_column not in frame_df.columns:
        raise ValueError(f"Missing required feature column: {args.feature_column}")

    triplets = _parse_triplets(args.triplets_txt)
    fold_map = _fold_map(args.susceptibility_xlsx)

    out_tables = args.output_dir / "tables"
    out_plots = args.output_dir / "plots"
    out_config = args.output_dir / "config"
    out_tables.mkdir(parents=True, exist_ok=True)
    out_plots.mkdir(parents=True, exist_ok=True)
    out_config.mkdir(parents=True, exist_ok=True)

    feature_slug = _feature_slug(args.feature_column)
    rep_df, mut_df = _mutation_summaries(frame_df, args.feature_column)
    rep_df.to_csv(out_tables / f"{feature_slug}_by_replicate.csv", index=False)
    mut_df.to_csv(out_tables / f"{feature_slug}_by_mutation.csv", index=False)

    triplet_summary_rows: list[dict[str, object]] = []
    for triplet in triplets:
        wanted = set(triplet)
        sub = frame_df[frame_df["mutation"].astype(str).isin(wanted)].copy()
        if sub["mutation"].nunique() != 3:
            raise ValueError(f"Missing mutations for triplet {triplet}")
        plot_name = f"{feature_slug}_triplet_{_safe_name(triplet[0])}_{_safe_name(triplet[1])}_{_safe_name(triplet[2])}.png"
        result = _plot_triplet(
            sub,
            rep_df,
            triplet,
            feature=args.feature_column,
            fold_map=fold_map,
            output_png=out_plots / plot_name,
        )
        triplet_summary_rows.extend(result["summary_rows"])

    pd.DataFrame(triplet_summary_rows).to_csv(out_tables / f"{feature_slug}_triplet_summary.csv", index=False)
    config = {
        "frame_feature_csv": str(args.frame_feature_csv),
        "triplets_txt": str(args.triplets_txt),
        "susceptibility_xlsx": str(args.susceptibility_xlsx),
        "feature_column": str(args.feature_column),
        "triplets": ["|".join(t) for t in triplets],
    }
    (out_config / "run_config.json").write_text(json.dumps(config, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
