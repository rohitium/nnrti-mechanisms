#!/usr/bin/env python3
"""Plot a standalone mean-occupancy heatmap across all simulated mutations."""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd

os.environ.setdefault("MPLBACKEND", "Agg")

from ..susceptibility import load_dor_susceptibilities
from .plot_dor_susceptibility_bars import CATEGORY_ORDER, NEGATIVE_CONTROLS, POSITIVE_CONTROLS, TEST_SET
from .plot_triplet_contact_story import (
    _aa1_to_aa3,
    _compute_triplet_contact_stats,
    _format_fold_label,
    _label_resname_for_position,
    _load_replicate_meta,
    _normalize_mutation_token,
    _parse_mutation_ops,
    _residue_label,
)


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plot a standalone all-mutation contact occupancy heatmap.")
    parser.add_argument("--manifest", type=Path, default=Path("manifests/md_manifest.csv"))
    parser.add_argument("--susceptibility-xlsx", type=Path, default=Path("data/DRM-susceptibilities.csv.xlsx"))
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results/analysis/triplet_story_analyses/contact_story_all_mutations_excluding_f227c"),
    )
    parser.add_argument("--ligand-resname", type=str, default="2KW")
    parser.add_argument("--resid-offset", type=int, default=-3)
    parser.add_argument("--contact-cutoff", type=float, default=4.0)
    parser.add_argument("--window-ns", type=float, default=100.0)
    parser.add_argument("--min-any-mean-occ-display", type=float, default=0.5)
    parser.add_argument("--exclude-mutations", type=str, default="F227C")
    parser.add_argument(
        "--output-prefix",
        type=str,
        default="all_mutation_mean_occupancy_heatmap_excluding_f227c",
    )
    return parser.parse_args()


def _fold_map(xlsx_path: Path) -> dict[str, float]:
    df = load_dor_susceptibilities(xlsx_path)
    out = {str(row["mutation"]): float(row["dor_fold_reduction"]) for _, row in df.iterrows()}
    out["WT"] = 1.0
    return out


def _mutation_order(manifest_csv: Path, fold_map: dict[str, float], excluded: set[str]) -> list[str]:
    mf = pd.read_csv(manifest_csv)
    muts = sorted({_normalize_mutation_token(m) for m in mf["mutation"].astype(str) if _normalize_mutation_token(m)})
    muts = [m for m in muts if m not in excluded]

    def _category_rank(mutation: str) -> int:
        if mutation == "WT":
            return -1
        if mutation in NEGATIVE_CONTROLS:
            return CATEGORY_ORDER["Negative control"]
        if mutation in POSITIVE_CONTROLS:
            return CATEGORY_ORDER["Positive control"]
        if mutation in TEST_SET:
            return CATEGORY_ORDER["Test set"]
        return 99

    return sorted(
        muts,
        key=lambda m: (_category_rank(m), float(fold_map.get(m, np.inf)), m),
    )


def _build_heatmap_matrix(
    mut_occ: pd.DataFrame,
    mutation_order: list[str],
    min_any_mean_occ_display: float,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    disp_filter = (
        mut_occ[mut_occ["mutation"].isin(mutation_order)]
        .groupby(["traj_resid", "auth_resid"], as_index=False)["occupancy_mean"]
        .max()
    )
    disp_filter["occupancy_mean"] = pd.to_numeric(disp_filter["occupancy_mean"], errors="coerce")
    disp_filter = disp_filter[disp_filter["occupancy_mean"] >= float(min_any_mean_occ_display)].copy()
    if disp_filter.empty:
        raise ValueError(
            f"No residues pass occupancy_mean >= {min_any_mean_occ_display:.2f} across the selected mutations."
        )

    residue_df = (
        disp_filter[["traj_resid", "auth_resid"]]
        .drop_duplicates()
        .sort_values(["auth_resid", "traj_resid"], kind="stable")
        .reset_index(drop=True)
    )

    def _wt_position_resname(auth_resid: int, traj_resid: int) -> str:
        pos = int(auth_resid)
        for mutation in mutation_order:
            for src, mpos, _dst in _parse_mutation_ops(mutation):
                if int(mpos) == pos:
                    return _aa1_to_aa3(src)
        return _label_resname_for_position(
            mut_occ=mut_occ,
            wt_mutation="WT",
            traj_resid=int(traj_resid),
            auth_resid=int(auth_resid),
        )

    residue_df["resname"] = residue_df.apply(
        lambda r: _wt_position_resname(auth_resid=int(r["auth_resid"]), traj_resid=int(r["traj_resid"])),
        axis=1,
    )
    residue_df["label"] = residue_df.apply(
        lambda r: _residue_label(int(r["auth_resid"]), str(r["resname"])),
        axis=1,
    )

    key_lookup = residue_df.copy()
    key_lookup["key"] = list(zip(key_lookup["traj_resid"], key_lookup["auth_resid"]))
    map_lab = dict(zip(key_lookup["key"], key_lookup["label"]))

    plot_df = mut_occ.copy()
    plot_df["key"] = list(zip(plot_df["traj_resid"], plot_df["auth_resid"]))
    plot_df = plot_df[plot_df["key"].isin(set(map_lab.keys()))].copy()
    plot_df["label"] = plot_df["key"].map(map_lab)

    hm = (
        plot_df.pivot_table(index="mutation", columns="label", values="occupancy_pooled", aggfunc="mean")
        .reindex(index=mutation_order, columns=residue_df["label"].tolist())
    )
    return residue_df, hm


def _plot_heatmap(matrix: pd.DataFrame, fold_map: dict[str, float], output_png: Path) -> None:
    import matplotlib.pyplot as plt
    from matplotlib.colors import BoundaryNorm

    arr = matrix.to_numpy(dtype=float)
    fig_w = min(22.0, 2.6 + 0.52 * arr.shape[1])
    fig_h = min(14.0, 1.8 + 0.42 * arr.shape[0])
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    cmap = plt.get_cmap("cividis")
    boundaries = np.arange(0.0, 1.0001, 0.025)
    norm = BoundaryNorm(boundaries=boundaries, ncolors=cmap.N, clip=True)
    im = ax.imshow(arr, aspect="auto", cmap=cmap, norm=norm)
    ax.set_yticks(
        np.arange(len(matrix.index)),
        [_format_fold_label(m, fold_map=fold_map) for m in matrix.index.tolist()],
    )
    ax.set_xticks(np.arange(len(matrix.columns)), matrix.columns.tolist(), rotation=45, ha="right")
    ax.set_xlabel("Residue")
    ax.set_ylabel("Mutation")
    ax.set_title("Mean Occupancy Heatmap")
    for i in range(arr.shape[0]):
        for j in range(arr.shape[1]):
            v = arr[i, j]
            if not np.isfinite(v):
                continue
            txt_color = "white" if float(v) < 0.62 else "black"
            ax.text(j, i, f"{float(v):.2f}", ha="center", va="center", fontsize=7, color=txt_color)
    cbar = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    cbar.set_ticks(np.arange(0.0, 1.01, 0.1))
    cbar.set_label("Mean occupancy")
    output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_png, dpi=320, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    args = _parse_args()
    out_plots = args.output_dir / "plots"
    out_tables = args.output_dir / "tables"
    out_config = args.output_dir / "config"
    out_plots.mkdir(parents=True, exist_ok=True)
    out_tables.mkdir(parents=True, exist_ok=True)
    out_config.mkdir(parents=True, exist_ok=True)

    excluded = {_normalize_mutation_token(x) for x in str(args.exclude_mutations).split(",") if _normalize_mutation_token(x)}
    fold_map = _fold_map(args.susceptibility_xlsx)
    mutation_order = _mutation_order(args.manifest, fold_map=fold_map, excluded=excluded)
    metas = _load_replicate_meta(args.manifest, needed_mutations=set(mutation_order))
    if not metas:
        raise ValueError("No replicate metadata found for selected mutations.")

    rep_occ, mut_occ, timing_df = _compute_triplet_contact_stats(
        metas=metas,
        mutation_triplet=tuple(mutation_order),
        ligand_resname=str(args.ligand_resname),
        resid_offset=int(args.resid_offset),
        contact_cutoff=float(args.contact_cutoff),
        window_ns=float(args.window_ns),
    )
    rep_occ.to_csv(out_tables / "replicate_contact.csv", index=False)
    mut_occ.to_csv(out_tables / "mutation_contact.csv", index=False)
    timing_df.to_csv(out_tables / "timing_audit.csv", index=False)

    residue_df, heatmap_df = _build_heatmap_matrix(
        mut_occ=mut_occ,
        mutation_order=mutation_order,
        min_any_mean_occ_display=float(args.min_any_mean_occ_display),
    )
    residue_df.to_csv(out_tables / "display_residues.csv", index=False)
    heatmap_df.reset_index().rename(columns={"index": "mutation"}).to_csv(
        out_tables / "mean_occupancy_heatmap_matrix.csv",
        index=False,
    )

    output_png = out_plots / f"{str(args.output_prefix)}.png"
    _plot_heatmap(heatmap_df, fold_map=fold_map, output_png=output_png)

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
                "min_any_mean_occ_display": float(args.min_any_mean_occ_display),
                "exclude_mutations": sorted(excluded),
                "n_mutations": int(len(mutation_order)),
                "n_display_residues": int(len(residue_df)),
                "output_png": str(output_png),
            },
            indent=2,
        )
    )
    print(f"Saved {output_png}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
