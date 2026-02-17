#!/usr/bin/env python3
"""Plot a key-contact occupancy heatmap for doravirine (DOR) across mutations.

Occupancy definition (default):
  distance_angstrom <= distance_ref_angstrom + margin_angstrom

Inputs:
  - results/dor_key_contacts_timeseries_by_mutation/*.csv
  - results/dor_key_contact_definitions_4ncg.csv
  - results/md_manifest.csv (for fold-reduction ordering)

Outputs:
  - results/plots/dor_key_contact_occupancy_heatmap.png (default)
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


_AA3_TO_1 = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
}


def _mutation_order(manifest_csv: Path) -> list[str]:
    m = pd.read_csv(manifest_csv)
    m["fold_reduction"] = pd.to_numeric(m.get("fold_reduction"), errors="coerce")
    uniq = (
        m[["mutation", "fold_reduction"]]
        .drop_duplicates()
        .sort_values(["fold_reduction", "mutation"], na_position="first")
    )
    muts = uniq["mutation"].astype(str).tolist()
    if "WT" in muts:
        muts = ["WT"] + [x for x in muts if x != "WT"]
    return muts


def _load_timeseries(timeseries_dir: Path) -> pd.DataFrame:
    paths = sorted(timeseries_dir.glob("*_dor_key_contacts_timeseries.csv"))
    if not paths:
        raise FileNotFoundError(f"No per-mutation key-contact timeseries CSVs found in {timeseries_dir}")
    parts: list[pd.DataFrame] = []
    for p in paths:
        parts.append(pd.read_csv(p))
    df = pd.concat(parts, ignore_index=True)
    return df


def _load_contact_defs(defs_csv: Path) -> pd.DataFrame:
    d = pd.read_csv(defs_csv).copy()
    d["contact_id"] = d["contact_id"].astype(str)
    d["category"] = d["category"].astype(str).str.lower()
    d["protein_resid_auth"] = pd.to_numeric(d["protein_resid_auth"], errors="coerce")
    d["protein_resname"] = d["protein_resname"].astype(str)
    d["protein_atom"] = d["protein_atom"].astype(str)
    d["ligand_atom"] = d["ligand_atom"].astype(str)
    d = d.sort_values(["category", "protein_resid_auth", "contact_id"], na_position="last")
    d["contact_short"] = d.apply(
        lambda r: f"{_AA3_TO_1.get(str(r['protein_resname']).upper(), str(r['protein_resname'])[:1])}"
        f"{int(r['protein_resid_auth']) if pd.notna(r['protein_resid_auth']) else '?'} "
        f"({r['protein_atom']}:{r['ligand_atom']})",
        axis=1,
    )
    return d


def _compute_occupancy_table(timeseries_dir: Path, margin_angstrom: float) -> pd.DataFrame:
    rep_occ = _compute_replicate_occupancy_table(timeseries_dir=timeseries_dir, margin_angstrom=margin_angstrom)
    occ = rep_occ.groupby("mutation", as_index=False).mean(numeric_only=True).set_index("mutation")
    return occ


def _compute_replicate_occupancy_table(timeseries_dir: Path, margin_angstrom: float) -> pd.DataFrame:
    df = _load_timeseries(timeseries_dir)
    if df.empty:
        raise ValueError("Empty key-contact timeseries table.")

    df["distance_angstrom"] = pd.to_numeric(df["distance_angstrom"], errors="coerce")
    df["distance_ref_angstrom"] = pd.to_numeric(df["distance_ref_angstrom"], errors="coerce")
    df = df.dropna(subset=["mutation", "replicate", "contact_id", "distance_angstrom", "distance_ref_angstrom"]).copy()

    thr = df["distance_ref_angstrom"] + float(margin_angstrom)
    df["is_contact"] = df["distance_angstrom"] <= thr

    rep_occ = (
        df.groupby(["mutation", "replicate", "contact_id"], as_index=False)["is_contact"]
        .mean()
        .rename(columns={"is_contact": "occupancy"})
    )
    rep_mat = (
        rep_occ.groupby(["mutation", "contact_id"], as_index=False)["occupancy"]
        .mean()
    )
    rep_wide = (
        rep_occ.pivot(index=["mutation", "replicate"], columns="contact_id", values="occupancy")
        .reset_index()
    )
    # Use stable mutation/contact coverage from grouped means to avoid sparse oddities.
    keep_contacts = rep_mat["contact_id"].astype(str).unique().tolist()
    cols = ["mutation", "replicate"] + [c for c in keep_contacts if c in rep_wide.columns]
    return rep_wide[cols]


def _compute_replicate_contact_deltas(timeseries_dir: Path, margin_angstrom: float, contact_ids: list[str]) -> pd.DataFrame:
    rep = _compute_replicate_occupancy_table(timeseries_dir=timeseries_dir, margin_angstrom=margin_angstrom).copy()
    contact_ids = [c for c in contact_ids if c in rep.columns]
    if not contact_ids:
        return pd.DataFrame()
    wt = rep[rep["mutation"].astype(str) == "WT"][["replicate"] + contact_ids].copy()
    if wt.empty:
        return pd.DataFrame()
    wt = wt.rename(columns={c: f"wt__{c}" for c in contact_ids})
    rows: list[dict[str, float | str]] = []
    for mut in sorted(rep["mutation"].astype(str).unique().tolist()):
        if mut == "WT":
            continue
        m = rep[rep["mutation"].astype(str) == mut][["replicate"] + contact_ids].copy()
        merged = m.merge(wt, on="replicate", how="inner")
        if merged.empty:
            continue
        row: dict[str, float | str] = {"mutation": mut, "n_reps": int(len(merged))}
        for cid in contact_ids:
            vals = pd.to_numeric(merged[cid], errors="coerce") - pd.to_numeric(merged[f"wt__{cid}"], errors="coerce")
            vals = vals.dropna().to_numpy(dtype=float)
            if vals.size == 0:
                row[f"{cid}__mean"] = np.nan
                row[f"{cid}__sem"] = np.nan
            else:
                row[f"{cid}__mean"] = float(np.nanmean(vals))
                row[f"{cid}__sem"] = float(np.nanstd(vals, ddof=1) / np.sqrt(vals.size)) if vals.size > 1 else 0.0
        rows.append(row)
    return pd.DataFrame(rows)


def plot_occupancy_heatmap(
    timeseries_dir: Path,
    contact_defs_csv: Path,
    manifest_csv: Path,
    output_png: Path,
    margin_angstrom: float,
    mode: str,
) -> None:
    import matplotlib.pyplot as plt
    from matplotlib.colors import TwoSlopeNorm

    occ = _compute_occupancy_table(timeseries_dir=timeseries_dir, margin_angstrom=margin_angstrom)
    defs = _load_contact_defs(contact_defs_csv)

    # Ordering.
    mut_order = _mutation_order(manifest_csv)
    mut_order = [m for m in mut_order if m in occ.index.astype(str)]
    contact_order = defs["contact_id"].astype(str).tolist()
    contact_order = [c for c in contact_order if c in occ.columns.astype(str)]
    occ = occ.reindex(index=mut_order, columns=contact_order)

    # Split contacts by category and render as disjoint panels.
    hydro = defs[defs["category"] == "hydrophobic"]["contact_id"].astype(str).tolist()
    polar = defs[defs["category"] == "polar"]["contact_id"].astype(str).tolist()
    hydro = [c for c in hydro if c in occ.columns]
    polar = [c for c in polar if c in occ.columns]
    plot_df = occ.copy()
    cmap = "viridis"
    norm = None
    title = "RT:DOR key contact ocupancy"
    cbar_label = "Occupancy"

    if mode == "delta":
        if "WT" not in occ.index.astype(str):
            raise ValueError("mode=delta requires WT row in occupancy table.")
        wt = occ.loc["WT"].to_numpy(dtype=float)
        wt_map = dict(zip(occ.columns.tolist(), wt.tolist()))
        plot_df = occ.copy()
        for cid in plot_df.columns:
            plot_df[cid] = pd.to_numeric(plot_df[cid], errors="coerce") - float(wt_map.get(cid, np.nan))
        cmap = "coolwarm"
        vals = plot_df.to_numpy(dtype=float)
        vmax = float(np.nanmax(np.abs(vals))) if np.isfinite(vals).any() else 1.0
        vmax = max(0.05, min(1.0, vmax))
        norm = TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax)
        title = "RT:DOR key contact ocupancy shift"
        cbar_label = "Δ occupancy (mutant - WT)"

    groups: list[tuple[str, list[str]]] = []
    if hydro:
        groups.append(("Hydrophobic", hydro))
    if polar:
        groups.append(("Polar", polar))
    if not groups:
        raise ValueError("No contacts available after ordering/filtering.")

    width_ratios = [max(1, len(cols)) for _, cols in groups]
    fig, axes = plt.subplots(
        1,
        len(groups),
        figsize=(min(18, 1.4 + 0.55 * sum(width_ratios)), 0.8 + 0.45 * len(mut_order)),
        sharey=True,
        gridspec_kw={"width_ratios": width_ratios, "wspace": 0.08},
    )
    axes = np.atleast_1d(axes)
    im_ref = None
    for i, (ax, (glabel, cols)) in enumerate(zip(axes, groups)):
        vals = plot_df.reindex(columns=cols).to_numpy(dtype=float)
        im = ax.imshow(vals, aspect="auto", interpolation="nearest", cmap=cmap, norm=norm)
        im_ref = im if im_ref is None else im_ref
        xtick_labels: list[str] = []
        for cid in cols:
            row = defs[defs["contact_id"] == cid].head(1)
            xtick_labels.append(str(row["contact_short"].iloc[0]) if not row.empty else cid)
        ax.set_xticks(range(len(cols)))
        ax.set_xticklabels(xtick_labels, rotation=45, ha="right", fontsize=8)
        ax.set_yticks(range(len(mut_order)))
        if i == 0:
            ax.set_yticklabels(mut_order, fontsize=9)
            ax.set_ylabel("Mutation")
        else:
            ax.tick_params(axis="y", which="both", left=False, labelleft=False)
        ax.set_xticks(np.arange(len(cols) + 1) - 0.5, minor=True)
        ax.set_yticks(np.arange(len(mut_order) + 1) - 0.5, minor=True)
        ax.grid(which="minor", color="white", linestyle="-", linewidth=1.0, alpha=0.8)
        ax.tick_params(which="minor", bottom=False, left=False)
        ax.text(
            0.5,
            -0.18,
            glabel,
            transform=ax.transAxes,
            ha="center",
            va="top",
            fontsize=10,
        )

    fig.suptitle(title, fontsize=12, fontweight="bold", y=0.98)
    cbar = fig.colorbar(im_ref, ax=axes.tolist(), fraction=0.026, pad=0.02)
    cbar.set_label(cbar_label)

    output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_png, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {output_png}")


def plot_selected_contact_correlations(
    timeseries_dir: Path,
    contact_defs_csv: Path,
    manifest_csv: Path,
    output_png: Path,
    margin_angstrom: float,
) -> None:
    import matplotlib.pyplot as plt
    from scipy import stats

    occ = _compute_occupancy_table(timeseries_dir=timeseries_dir, margin_angstrom=margin_angstrom)
    manifest = pd.read_csv(manifest_csv)[["mutation", "fold_reduction"]].drop_duplicates().copy()
    manifest["fold_reduction"] = pd.to_numeric(manifest["fold_reduction"], errors="coerce")

    targets = [
        ("polar_LYS103_N_N19", "K103 (N:N19)"),
        ("hydrophobic_TYR181_CD1_F14", "Y181 (CD1:F14)"),
        ("hydrophobic_TYR188_CB_F14", "Y188 (CB:F14)"),
    ]
    available_targets = [(cid, label) for cid, label in targets if cid in occ.columns]
    if not available_targets:
        raise ValueError("None of requested target contacts were found in occupancy table.")

    contact_ids = [cid for cid, _ in available_targets]
    d = _compute_replicate_contact_deltas(
        timeseries_dir=timeseries_dir,
        margin_angstrom=margin_angstrom,
        contact_ids=contact_ids,
    )
    if d.empty:
        raise ValueError("Could not compute replicate-matched occupancy deltas.")
    d = d.merge(manifest.rename(columns={"fold_reduction": "fold_reduction"}), on="mutation", how="left")
    d = d[pd.to_numeric(d["fold_reduction"], errors="coerce") > 0].copy()
    d["fold_reduction"] = pd.to_numeric(d["fold_reduction"], errors="coerce")
    d["log10_fc"] = np.log10(d["fold_reduction"])
    if d.empty:
        raise ValueError("No mutant rows with valid fold-reduction found.")
    d["is_combo"] = d["mutation"].astype(str).str.contains(r"\+")

    fig, axes = plt.subplots(1, len(available_targets), figsize=(5.6 * len(available_targets), 4.6), constrained_layout=True)
    axes = np.atleast_1d(axes)
    for ax, (cid, short_label) in zip(axes, available_targets):
        singles = d[~d["is_combo"]]
        combos = d[d["is_combo"]]
        y_col = f"{cid}__mean"
        e_col = f"{cid}__sem"
        ax.errorbar(
            singles["log10_fc"],
            singles[y_col],
            yerr=singles[e_col] if e_col in singles.columns else None,
            fmt="o",
            color="#1f77b4",
            markersize=6,
            capsize=3,
            alpha=0.9,
            label="Single DRM",
        )
        ax.errorbar(
            combos["log10_fc"],
            combos[y_col],
            yerr=combos[e_col] if e_col in combos.columns else None,
            fmt="s",
            color="#d62728",
            markersize=5.5,
            capsize=3,
            alpha=0.9,
            label="DRM Combination",
        )

        for _, r in d.iterrows():
            ax.annotate(
                r["mutation"],
                (r["log10_fc"], r[y_col]),
                textcoords="offset points",
                xytext=(5, 4),
                fontsize=7,
                bbox=dict(boxstyle="round,pad=0.2", facecolor="white", alpha=0.7, edgecolor="none"),
            )

        ax.axhline(0.0, color="black", linewidth=1.1, alpha=0.9, label="WT")
        x = d["log10_fc"].to_numpy(dtype=float)
        y = d[y_col].to_numpy(dtype=float)
        if np.isfinite(x).sum() >= 3 and np.isfinite(y).sum() >= 3 and np.nanstd(x) > 1e-12 and np.nanstd(y) > 1e-12:
            slope, intercept = np.polyfit(x, y, 1)
            xline = np.linspace(np.nanmin(x), np.nanmax(x), 100)
            yline = slope * xline + intercept
            ax.plot(xline, yline, "--", color="#777777", linewidth=1.7, alpha=0.9)
            r, p = stats.pearsonr(x, y)
            trend_lbl = f"Trend line (R^2={r*r:.3f}, p={p:.3g})"
        else:
            trend_lbl = "Trend line"
        ax.plot([], [], "--", color="#777777", linewidth=1.7, label=trend_lbl)

        ax.set_title(f"{short_label} vs " + r"$\log_{10}(\mathrm{Fold\ Change})$", fontsize=11, fontweight="normal")
        ax.set_xlabel(r"$\log_{10}(\mathrm{Fold\ Change})$")
        ax.set_ylabel("Δ occupancy (mutant - WT)")
        ax.grid(alpha=0.3, linestyle=":", linewidth=0.8)
        ax.legend(loc="upper left", frameon=True, fontsize=8)

    output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_png, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {output_png}")


def plot_all_contact_correlations(
    timeseries_dir: Path,
    contact_defs_csv: Path,
    manifest_csv: Path,
    output_png: Path,
    margin_angstrom: float,
) -> None:
    import matplotlib.pyplot as plt
    from scipy import stats

    occ = _compute_occupancy_table(timeseries_dir=timeseries_dir, margin_angstrom=margin_angstrom)
    defs = _load_contact_defs(contact_defs_csv)
    manifest = pd.read_csv(manifest_csv)[["mutation", "fold_reduction"]].drop_duplicates().copy()
    manifest["fold_reduction"] = pd.to_numeric(manifest["fold_reduction"], errors="coerce")

    contact_ids = defs["contact_id"].astype(str).tolist()
    contact_ids = [cid for cid in contact_ids if cid in occ.columns]
    if not contact_ids:
        raise ValueError("No contact IDs found for all-contact correlation plot.")

    d = _compute_replicate_contact_deltas(
        timeseries_dir=timeseries_dir,
        margin_angstrom=margin_angstrom,
        contact_ids=contact_ids,
    )
    if d.empty:
        raise ValueError("Could not compute replicate-matched occupancy deltas for all contacts.")
    d = d.merge(manifest, on="mutation", how="left")
    d = d[pd.to_numeric(d["fold_reduction"], errors="coerce") > 0].copy()
    d["fold_reduction"] = pd.to_numeric(d["fold_reduction"], errors="coerce")
    d["log10_fc"] = np.log10(d["fold_reduction"])
    d["is_combo"] = d["mutation"].astype(str).str.contains(r"\+")

    n = len(contact_ids)
    ncols = 3
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(5.1 * ncols, 3.8 * nrows), constrained_layout=True)
    axes = np.atleast_1d(axes).flatten()

    for i, cid in enumerate(contact_ids):
        ax = axes[i]
        label_row = defs[defs["contact_id"] == cid].head(1)
        short_label = str(label_row["contact_short"].iloc[0]) if not label_row.empty else cid

        y_col = f"{cid}__mean"
        e_col = f"{cid}__sem"
        if y_col not in d.columns:
            ax.set_visible(False)
            continue

        singles = d[~d["is_combo"]]
        combos = d[d["is_combo"]]
        ax.errorbar(
            singles["log10_fc"],
            singles[y_col],
            yerr=singles[e_col] if e_col in singles.columns else None,
            fmt="o",
            color="#1f77b4",
            markersize=5,
            capsize=2.5,
            alpha=0.9,
            label="Single DRM",
        )
        ax.errorbar(
            combos["log10_fc"],
            combos[y_col],
            yerr=combos[e_col] if e_col in combos.columns else None,
            fmt="s",
            color="#d62728",
            markersize=4.6,
            capsize=2.5,
            alpha=0.9,
            label="DRM Combination",
        )

        ax.axhline(0.0, color="black", linewidth=1.0, alpha=0.9, label="WT")
        x = d["log10_fc"].to_numpy(dtype=float)
        y = d[y_col].to_numpy(dtype=float)
        if np.isfinite(x).sum() >= 3 and np.isfinite(y).sum() >= 3 and np.nanstd(x) > 1e-12 and np.nanstd(y) > 1e-12:
            slope, intercept = np.polyfit(x, y, 1)
            xline = np.linspace(np.nanmin(x), np.nanmax(x), 100)
            yline = slope * xline + intercept
            ax.plot(xline, yline, "--", color="#777777", linewidth=1.4, alpha=0.9)
            r, p = stats.pearsonr(x, y)
            trend_lbl = f"Trend line (R^2={r*r:.3f}, p={p:.3g})"
        else:
            trend_lbl = "Trend line"
        ax.plot([], [], "--", color="#777777", linewidth=1.4, label=trend_lbl)

        ax.set_title(short_label, fontsize=10, fontweight="normal")
        ax.set_xlabel(r"$\log_{10}(\mathrm{Fold\ Change})$", fontsize=9)
        ax.set_ylabel("Δ occupancy (mutant - WT)", fontsize=9)
        ax.grid(alpha=0.3, linestyle=":", linewidth=0.8)
        ax.legend(loc="upper left", frameon=True, fontsize=7)

    for j in range(n, len(axes)):
        axes[j].set_visible(False)

    output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_png, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {output_png}")


def main() -> int:
    parser = argparse.ArgumentParser(description="Plot doravirine key-contact occupancy heatmap.")
    parser.add_argument("--timeseries-dir", type=Path, default=Path("results/dor_key_contacts_timeseries_by_mutation"))
    parser.add_argument("--contact-defs", type=Path, default=Path("results/dor_key_contact_definitions_4ncg.csv"))
    parser.add_argument("--manifest", type=Path, default=Path("results/md_manifest.csv"))
    parser.add_argument("--output", type=Path, default=Path("results/plots/dor_key_contact_occupancy_heatmap.png"))
    parser.add_argument("--margin-angstrom", type=float, default=1.0)
    parser.add_argument("--mode", choices=["absolute", "delta"], default="delta")
    parser.add_argument(
        "--corr-output",
        type=Path,
        default=Path("results/plots/dor_key_contact_selected_vs_fold_reduction.png"),
    )
    parser.add_argument(
        "--corr-all-output",
        type=Path,
        default=Path("results/plots/dor_key_contact_all_vs_fold_reduction.png"),
    )
    args = parser.parse_args()

    for p in [args.timeseries_dir, args.contact_defs, args.manifest]:
        if not p.exists():
            raise FileNotFoundError(f"Missing: {p}")

    plot_occupancy_heatmap(
        timeseries_dir=args.timeseries_dir,
        contact_defs_csv=args.contact_defs,
        manifest_csv=args.manifest,
        output_png=args.output,
        margin_angstrom=float(args.margin_angstrom),
        mode=str(args.mode),
    )
    plot_selected_contact_correlations(
        timeseries_dir=args.timeseries_dir,
        contact_defs_csv=args.contact_defs,
        manifest_csv=args.manifest,
        output_png=args.corr_output,
        margin_angstrom=float(args.margin_angstrom),
    )
    plot_all_contact_correlations(
        timeseries_dir=args.timeseries_dir,
        contact_defs_csv=args.contact_defs,
        manifest_csv=args.manifest,
        output_png=args.corr_all_output,
        margin_angstrom=float(args.margin_angstrom),
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
