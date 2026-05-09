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
from .plot_dor_susceptibility_bars import CATEGORY_COLORS, CATEGORY_ORDER, NEGATIVE_CONTROLS, POSITIVE_CONTROLS, TEST_SET
from .plot_triplet_contact_story import (
    _aa1_to_aa3,
    _compute_triplet_contact_stats,
    _format_fold_label,
    _label_resname_for_position,
    _load_replicate_meta,
    _normalize_mutation_token,
    _parse_mutation_ops,
    _residue_label,
    _wt_contact_region,
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
    parser.add_argument(
        "--include-all-contacted-residues",
        action="store_true",
        help="Display every residue that contacted DOR at least once in any selected mutation.",
    )
    parser.add_argument(
        "--wt-reference",
        action="store_true",
        help="Plot occupancy differences relative to WT for each residue.",
    )
    parser.add_argument(
        "--drop-wt-row",
        action="store_true",
        help="Drop the WT row from a WT-referenced heatmap.",
    )
    parser.add_argument(
        "--restrict-to-wt-contacted-residues",
        action="store_true",
        help="Restrict displayed residues to the WT all-contacted-residue panel.",
    )
    parser.add_argument(
        "--wt-contact-table",
        type=Path,
        default=Path("results/analysis/triplet_story_analyses/contact_story_100ns/tables/mutation_contact.csv"),
        help="WT contact table used to define the original all-contacted-residue panel.",
    )
    parser.add_argument(
        "--group-by-wt-contact-region",
        action="store_true",
        help="Group displayed residues using the 7 WT NNIBP contact regions.",
    )
    parser.add_argument(
        "--group-by-mutation-set",
        action="store_true",
        help="Group heatmap rows by negative control, positive control, and limited-data mutation sets.",
    )
    parser.add_argument(
        "--reuse-existing-tables",
        action="store_true",
        help="Reuse mutation_contact.csv/replicate_contact.csv in the output directory instead of recomputing contacts.",
    )
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


def _load_wt_contact_keys(wt_contact_table: Path) -> set[tuple[int, int]]:
    df = pd.read_csv(wt_contact_table)
    missing = {"mutation", "traj_resid", "auth_resid"}.difference(df.columns)
    if missing:
        raise ValueError(f"{wt_contact_table} is missing required columns: {sorted(missing)}")
    wt_df = df[df["mutation"].astype(str).map(_normalize_mutation_token).eq("WT")].copy()
    if wt_df.empty:
        raise ValueError(f"No WT rows found in {wt_contact_table}.")
    wt_df["traj_resid"] = pd.to_numeric(wt_df["traj_resid"], errors="coerce")
    wt_df["auth_resid"] = pd.to_numeric(wt_df["auth_resid"], errors="coerce")
    wt_df = wt_df.dropna(subset=["traj_resid", "auth_resid"])
    return {(int(r.traj_resid), int(r.auth_resid)) for r in wt_df.itertuples(index=False)}


def _mutation_set_label(mutation: str) -> str:
    mutation = _normalize_mutation_token(mutation)
    if mutation in NEGATIVE_CONTROLS:
        return "Negative control"
    if mutation in POSITIVE_CONTROLS:
        return "Positive control"
    if mutation in TEST_SET:
        return "Limited data"
    if mutation == "WT":
        return "WT"
    return "Other"


def _mutation_set_color(label: str) -> str:
    if label == "Limited data":
        return CATEGORY_COLORS["Test set"]
    return CATEGORY_COLORS.get(label, "#9d9da1")


def _build_heatmap_matrix(
    mut_occ: pd.DataFrame,
    mutation_order: list[str],
    min_any_mean_occ_display: float,
    include_all_contacted_residues: bool,
    wt_reference: bool,
    drop_wt_row: bool,
    restrict_keys: set[tuple[int, int]] | None,
    group_by_wt_contact_region: bool,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    disp_filter = (
        mut_occ[mut_occ["mutation"].isin(mutation_order)]
        .groupby(["traj_resid", "auth_resid"], as_index=False)["occupancy_mean"]
        .max()
    )
    disp_filter["occupancy_mean"] = pd.to_numeric(disp_filter["occupancy_mean"], errors="coerce")
    if restrict_keys is not None:
        disp_filter["key"] = list(zip(disp_filter["traj_resid"], disp_filter["auth_resid"]))
        disp_filter = disp_filter[disp_filter["key"].isin(restrict_keys)].drop(columns=["key"]).copy()
    elif not include_all_contacted_residues:
        disp_filter = disp_filter[disp_filter["occupancy_mean"] >= float(min_any_mean_occ_display)].copy()
    if disp_filter.empty:
        if restrict_keys is not None:
            raise ValueError("None of the WT contacted residues were found across the selected mutations.")
        if include_all_contacted_residues:
            raise ValueError("No contacted residues were found across the selected mutations.")
        raise ValueError(f"No residues pass occupancy_mean >= {min_any_mean_occ_display:.2f} across the selected mutations.")

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
    if group_by_wt_contact_region:
        region_info = residue_df["auth_resid"].map(lambda x: _wt_contact_region(int(x)))
        residue_df["region_order"] = region_info.map(lambda x: int(x[0]))
        residue_df["pocket_region"] = region_info.map(lambda x: str(x[1]))
        residue_df["region_color"] = region_info.map(lambda x: str(x[2]))
        residue_df = residue_df.sort_values(
            ["region_order", "auth_resid", "traj_resid"],
            kind="stable",
        ).reset_index(drop=True)

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
    if include_all_contacted_residues or wt_reference:
        hm = hm.fillna(0.0)
    if wt_reference:
        if "WT" not in hm.index:
            raise ValueError("WT row is required for --wt-reference.")
        wt = hm.loc["WT"].fillna(0.0)
        hm = hm.subtract(wt, axis="columns")
    if drop_wt_row:
        hm = hm.drop(index="WT", errors="ignore")
    return residue_df, hm


def _plot_heatmap(
    matrix: pd.DataFrame,
    residue_df: pd.DataFrame,
    fold_map: dict[str, float],
    output_png: Path,
    wt_reference: bool,
    group_by_wt_contact_region: bool,
    group_by_mutation_set: bool,
) -> None:
    import matplotlib.pyplot as plt
    from matplotlib.colors import BoundaryNorm
    from matplotlib.transforms import blended_transform_factory
    import textwrap

    arr = matrix.to_numpy(dtype=float)
    fig_w = min(40.0, 2.8 + 0.78 * arr.shape[1])
    fig_h = min(14.0, 1.8 + 0.42 * arr.shape[0])
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    if wt_reference:
        cmap = plt.get_cmap("coolwarm")
        boundaries = np.arange(-1.0, 1.0001, 0.05)
    else:
        cmap = plt.get_cmap("cividis")
        boundaries = np.arange(0.0, 1.0001, 0.025)
    norm = BoundaryNorm(boundaries=boundaries, ncolors=cmap.N, clip=True)
    im = ax.imshow(arr, aspect="auto", cmap=cmap, norm=norm)
    ax.set_yticks(
        np.arange(len(matrix.index)),
        [_format_fold_label(m, fold_map=fold_map) for m in matrix.index.tolist()],
    )
    ax.set_xticks(np.arange(len(matrix.columns)), matrix.columns.tolist(), rotation=45, ha="right")
    ax.tick_params(axis="x", labelsize=15)
    ax.tick_params(axis="y", labelsize=15)
    ax.set_xlabel("Residue", fontsize=17)
    ax.set_ylabel("Mutation", fontsize=17)
    ax.set_title(
        "WT-referenced DOR Contact Occupancy Heatmap" if wt_reference else "Mean Occupancy Heatmap",
        fontsize=19,
        pad=84 if group_by_wt_contact_region else 12,
    )
    if group_by_wt_contact_region and {"pocket_region", "region_color"}.issubset(residue_df.columns):
        trans = blended_transform_factory(ax.transData, ax.transAxes)
        for _region, g in residue_df.groupby("pocket_region", sort=False):
            left = int(g.index.min())
            right = int(g.index.max())
            center = (left + right) / 2.0
            color = str(g["region_color"].iloc[0])
            label = "\n".join(textwrap.wrap(str(_region), width=16))
            ax.plot([left - 0.45, right + 0.45], [1.015, 1.015], color=color, linewidth=5, transform=trans, clip_on=False)
            ax.text(center, 1.045, label, ha="center", va="bottom", fontsize=15, fontweight="bold", color=color, transform=trans, clip_on=False)
            if right < len(residue_df) - 1:
                ax.axvline(right + 0.5, color="white", linewidth=2.0)
    if group_by_mutation_set:
        row_sets = pd.Series([_mutation_set_label(m) for m in matrix.index], index=matrix.index)
        ytrans = blended_transform_factory(ax.transAxes, ax.transData)
        for set_label, rows in row_sets.groupby(row_sets, sort=False):
            row_positions = [matrix.index.get_loc(idx) for idx in rows.index]
            top = min(row_positions) - 0.5
            bottom = max(row_positions) + 0.5
            center = (top + bottom) / 2.0
            color = _mutation_set_color(str(set_label))
            ax.plot([1.012, 1.012], [top, bottom], color=color, linewidth=8, solid_capstyle="butt", transform=ytrans, clip_on=False)
            ax.text(
                1.04,
                center,
                str(set_label).replace(" ", "\n"),
                ha="center",
                va="center",
                rotation=90,
                fontsize=15,
                fontweight="bold",
                color=color,
                transform=ytrans,
                clip_on=False,
            )
            if bottom < len(matrix.index) - 0.5:
                ax.axhline(bottom, color="white", linewidth=2.0)
    for i in range(arr.shape[0]):
        for j in range(arr.shape[1]):
            v = arr[i, j]
            if not np.isfinite(v):
                continue
            if wt_reference:
                txt_color = "white" if abs(float(v)) > 0.45 else "black"
            else:
                txt_color = "white" if float(v) < 0.62 else "black"
            ax.text(j, i, f"{float(v):.2f}", ha="center", va="center", fontsize=10, color=txt_color)
    cbar = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.08 if group_by_mutation_set else 0.02)
    if wt_reference:
        cbar.set_ticks(np.arange(-1.0, 1.01, 0.2))
        cbar.set_label("Δ occupancy vs WT", fontsize=15)
    else:
        cbar.set_ticks(np.arange(0.0, 1.01, 0.1))
        cbar.set_label("Mean occupancy", fontsize=15)
    cbar.ax.tick_params(labelsize=14)
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
    restrict_keys = _load_wt_contact_keys(args.wt_contact_table) if args.restrict_to_wt_contacted_residues else None
    if args.reuse_existing_tables:
        mut_occ_path = out_tables / "mutation_contact.csv"
        rep_occ_path = out_tables / "replicate_contact.csv"
        if not mut_occ_path.exists():
            raise FileNotFoundError(f"--reuse-existing-tables requested, but {mut_occ_path} does not exist.")
        mut_occ = pd.read_csv(mut_occ_path)
        rep_occ = pd.read_csv(rep_occ_path) if rep_occ_path.exists() else pd.DataFrame()
        timing_df = pd.DataFrame()
    else:
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
        include_all_contacted_residues=bool(args.include_all_contacted_residues),
        wt_reference=bool(args.wt_reference),
        drop_wt_row=bool(args.drop_wt_row),
        restrict_keys=restrict_keys,
        group_by_wt_contact_region=bool(args.group_by_wt_contact_region),
    )
    default_prefix = "all_mutation_mean_occupancy_heatmap_excluding_f227c"
    if str(args.output_prefix) == default_prefix:
        residue_table = out_tables / "display_residues.csv"
        matrix_table = out_tables / "mean_occupancy_heatmap_matrix.csv"
    else:
        residue_table = out_tables / f"{str(args.output_prefix)}_display_residues.csv"
        matrix_table = out_tables / f"{str(args.output_prefix)}_matrix.csv"
    residue_df.to_csv(residue_table, index=False)
    heatmap_df.reset_index().rename(columns={"index": "mutation"}).to_csv(matrix_table, index=False)

    output_png = out_plots / f"{str(args.output_prefix)}.png"
    _plot_heatmap(
        heatmap_df,
        residue_df=residue_df,
        fold_map=fold_map,
        output_png=output_png,
        wt_reference=bool(args.wt_reference),
        group_by_wt_contact_region=bool(args.group_by_wt_contact_region),
        group_by_mutation_set=bool(args.group_by_mutation_set),
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
                "min_any_mean_occ_display": float(args.min_any_mean_occ_display),
                "include_all_contacted_residues": bool(args.include_all_contacted_residues),
                "wt_reference": bool(args.wt_reference),
                "drop_wt_row": bool(args.drop_wt_row),
                "restrict_to_wt_contacted_residues": bool(args.restrict_to_wt_contacted_residues),
                "wt_contact_table": str(args.wt_contact_table),
                "group_by_wt_contact_region": bool(args.group_by_wt_contact_region),
                "group_by_mutation_set": bool(args.group_by_mutation_set),
                "reuse_existing_tables": bool(args.reuse_existing_tables),
                "exclude_mutations": sorted(excluded),
                "n_mutations": int(len(mutation_order)),
                "n_display_residues": int(len(residue_df)),
                "output_png": str(output_png),
                "residue_table": str(residue_table),
                "matrix_table": str(matrix_table),
            },
            indent=2,
        )
    )
    print(f"Saved {output_png}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
