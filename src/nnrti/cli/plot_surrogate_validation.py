#!/usr/bin/env python3
"""Do the intra-protein surrogates reproduce the Table 3 contact metrics?

Seven of the nine Table 3 columns measure a protein-ligand contact, so they
vanish in the unliganded enzyme and the apo/holo comparison cannot be run on
them. Each is therefore replaced by an intra-protein distance that spans the
pocket along the same direction, which exists with or without doravirine.

This figure is the evidence that the replacement is faithful: per genotype, the
surrogate against the Table 3 metric it stands in for, in the holo state where
both are defined. One panel per surrogate, drawn in the style of Figure 2D.

Every anchor atom is a backbone atom or a Calpha where the anchor residue is
mutated somewhere in the panel, so no surrogate can move merely because a side
chain was replaced.

Usage
-----
    PYTHONPATH=src python -m nnrti.cli.plot_surrogate_validation
"""
from __future__ import annotations

import argparse
import logging
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from nnrti.analysis.panel import RESISTANT, SUSCEPTIBLE, UNCERTAIN
from nnrti.fep.combine_neq import _repel_labels

REPO_ROOT = Path(__file__).resolve().parents[3]

COLORS = {"Susceptible": "#2E7D32", "Resistant": "#C62828", "Uncertain": "#9E9E9E"}
MARKERS = {"Susceptible": "o", "Resistant": "s", "Uncertain": "^"}
WT_COLOR = "#1a1a1a"

AXIS_LABEL_FONTSIZE = 15
TICK_LABEL_FONTSIZE = 12
POINT_LABEL_FONTSIZE = 8
LEGEND_FONTSIZE = 11
STAT_FONTSIZE = 11

# How each Table 3 column reads on an axis, and what the surrogate measures.
PRETTY = {
    "d(103 C=O...N4x)": r"$d$(Lys103 C=O$\cdots$N4x)",
    "d(S105...DOR)":    r"$d$(Ser105$\cdots$DOR)",
    "d(106...DOR)":     r"$d$(residue 106$\cdots$DOR)",
    "d(V179...DOR)":    r"$d$(Val179$\cdots$DOR)",
    "d(103 sc...DOR)":  r"$d$(Lys103 side chain$\cdots$DOR)",
}
ANCHOR_LABEL = {
    "d(103 C=O...N4x)": ("Lys103 N/O", "as_table3"),
    "d(103 sc...DOR)":  (r"Lys103 C$\alpha$", "ca"),
    "d(106...DOR)":     ("residue 106 side chain", "sc"),
    "d(S105...DOR)":    ("Ser105", "as_table3"),
    "d(V179...DOR)":    ("Val179", "as_table3"),
}
THREE = {"ALA":"Ala","ARG":"Arg","ASN":"Asn","ASP":"Asp","CYS":"Cys","GLN":"Gln",
         "GLU":"Glu","GLY":"Gly","HIS":"His","ILE":"Ile","LEU":"Leu","LYS":"Lys",
         "MET":"Met","PHE":"Phe","PRO":"Pro","SER":"Ser","THR":"Thr","TRP":"Trp",
         "TYR":"Tyr","VAL":"Val"}


def surrogate_label(column: str, partner: str) -> str:
    """'PRO226:C' + its column -> 'd(Lys103 N/O...Pro226 C)'."""
    anchor = ANCHOR_LABEL.get(column, (column, ""))[0]
    resname, _, atom = partner.partition(":")
    num = "".join(c for c in resname if c.isdigit())
    name = THREE.get("".join(c for c in resname if c.isalpha()).upper(), resname)
    return rf"$d$({anchor}$\cdots${name}{num} {atom})"


def category_for(genotype: str) -> str:
    if genotype == "WT":
        return "Wild type"
    if genotype in SUSCEPTIBLE:
        return "Susceptible"
    if genotype in RESISTANT:
        return "Resistant"
    return "Uncertain"


def main(argv: list[str] | None = None) -> int:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--surrogate-frames", type=Path,
                    default=REPO_ROOT / "results/analysis/pocket_surrogates/tables/"
                                        "surrogate_per_frame_both.csv")
    ap.add_argument("--mechanisms", type=Path,
                    default=REPO_ROOT / "results/analysis/mechanisms/"
                                        "mechanism_coordinates.csv")
    ap.add_argument("--selected", type=Path,
                    default=REPO_ROOT / "results/analysis/pocket_surrogates/tables/"
                                        "surrogates_selected.csv")
    ap.add_argument("--output", type=Path,
                    default=REPO_ROOT / "results/analysis/pocket_surrogates/"
                                        "surrogate_validation.png")
    ap.add_argument("--table3", type=Path,
                    default=REPO_ROOT / "paper/tables/Table-3-structural.csv",
                    help="Published Table 3; its V(NNIBP) column is used directly.")
    ap.add_argument("--dpi", type=int, default=300)
    args = ap.parse_args(argv)

    sel = pd.read_csv(args.selected)

    # V(NNIBP) is a property of the protein alone, so it already exists in the
    # unliganded state and needs no stand-in. It is drawn here so the figure
    # shows every feature the apo/holo comparison uses; both axes are the same
    # quantity by construction, which is the point.
    direct = sel[sel.surrogate == "DIRECT"] if "surrogate" in sel else sel.iloc[0:0]
    sel = sel[sel.surrogate != "DIRECT"] if "surrogate" in sel else sel
    sur = pd.read_csv(args.surrogate_frames)
    sur = sur[sur.state == "holo"]          # both metrics are defined only here
    mech = pd.read_csv(args.mechanisms)

    coord_of = dict(zip(sel.table3_column, sel.table3_column))
    from nnrti.cli.find_pocket_surrogates import AXES
    coord_of = {c: co for c, co, _, _, _ in AXES}

    keys = ["mutation", "replicate", "frame"]
    coords = sorted({coord_of[c] for c in sel.table3_column})
    merged = sur.merge(mech[keys + coords], on=keys, how="inner")
    # replicate means first, then the genotype mean -- the way Table 3 is built
    per_rep = merged.groupby(["mutation", "replicate"]).mean(numeric_only=True).reset_index()
    gm = per_rep.groupby("mutation").mean(numeric_only=True)
    sem = per_rep.groupby("mutation").sem(numeric_only=True)

    vol = None
    if len(direct):
        # Table 3's own V(NNIBP) column, as published ("228.4 ± 11.7"). The
        # per-replicate pocket_volume_proxy in ddg_full.csv is computed over an
        # inconsistent number of frames per replicate (28 to 100) and puts
        # F227C at 540 +/- 300 A^3, which is not a Table 3 genotype at all.
        t3 = pd.read_csv(args.table3)
        col = next(c for c in t3.columns if "V(NNIBP)" in c)
        parsed = t3[col].astype(str).str.split("±", expand=True)
        # .to_numpy(), not the Series: passing a string index to the DataFrame
        # constructor makes pandas reindex by those labels, giving all NaN.
        vol = pd.DataFrame({
            "mean": pd.to_numeric(parsed[0].str.strip(), errors="coerce").to_numpy(),
            "sem": pd.to_numeric(parsed[1].str.strip(), errors="coerce").to_numpy(),
        }, index=t3["Genotype"].astype(str).to_numpy()).dropna(subset=["mean"])

    n = len(sel) + (1 if vol is not None else 0)
    ncol = 2 if n <= 4 else 3
    nrow = int(np.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(7.0 * ncol, 6.0 * nrow), squeeze=False)
    stats_rows = []

    for ax, (_, row) in zip(axes.ravel(), sel.iterrows()):
        column, feature = row.table3_column, row.surrogate
        coord = coord_of[column]
        sub = gm[[feature, coord]].dropna()
        xs, ys = sub[feature].to_numpy(), sub[coord].to_numpy()
        names = list(sub.index)

        lr = stats.linregress(xs, ys)
        rho, rho_p = stats.spearmanr(xs, ys)
        # Leave-one-out on R2, so a panel carried by a single high-leverage
        # genotype is visible in the table rather than only in the figure.
        loo = []
        for i in range(len(xs)):
            k = np.ones(len(xs), bool); k[i] = False
            loo.append((stats.linregress(xs[k], ys[k]).rvalue ** 2, names[i]))
        loo_min, loo_min_g = min(loo)
        stats_rows.append(dict(table3_column=column, surrogate=feature, n=len(xs),
                               r2=round(lr.rvalue ** 2, 3), p_value=lr.pvalue,
                               r2_min_leave_one_out=round(loo_min, 3),
                               most_influential=loo_min_g,
                               spearman_rho=round(rho, 3), spearman_p=rho_p))

        grid = np.linspace(xs.min(), xs.max(), 50)
        ax.plot(grid, lr.intercept + lr.slope * grid, "-", color="0.35",
                lw=1.4, zorder=2)

        seen = set()
        for g, x, y in zip(names, xs, ys):
            cat = category_for(g)
            ax.errorbar(x, y,
                        xerr=sem.loc[g, feature] if g in sem.index else None,
                        yerr=sem.loc[g, coord] if g in sem.index else None,
                        fmt=MARKERS.get(cat, "o"),
                        color=WT_COLOR if cat == "Wild type" else COLORS[cat],
                        ecolor="0.7", elinewidth=0.9, capsize=2, ms=10 if cat == "Wild type" else 7,
                        zorder=4 if cat == "Wild type" else 3,
                        label=cat if cat not in seen else None)
            seen.add(cat)

        _repel_labels(ax, fig, list(xs), list(ys), names,
                      fontsize=POINT_LABEL_FONTSIZE, n_iter=1200)
        ax.margins(0.13)

        ax.set_xlabel(f"{surrogate_label(column, row.partner)}  (Å)",
                      fontsize=AXIS_LABEL_FONTSIZE)
        ax.set_ylabel(f"{PRETTY.get(column, column)}  (Å)",
                      fontsize=AXIS_LABEL_FONTSIZE)
        ax.tick_params(axis="both", labelsize=TICK_LABEL_FONTSIZE)
        pstr = "$P$ < 0.001" if lr.pvalue < 1e-3 else f"$P$ = {lr.pvalue:.3f}"
        ax.annotate(f"$R^2$ = {lr.rvalue**2:.2f},  {pstr},  $n$ = {len(xs)}",
                    xy=(0.5, 1.015), xycoords="axes fraction", va="bottom", ha="center",
                    fontsize=STAT_FONTSIZE)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)

    if vol is not None:
        ax = axes.ravel()[len(sel)]
        xs = ys = vol["mean"].to_numpy()
        names = list(vol.index)
        lo, hi = xs.min(), xs.max()
        pad = 0.06 * (hi - lo)
        ax.plot([lo - pad, hi + pad], [lo - pad, hi + pad], "-", color="0.35",
                lw=1.4, zorder=2)
        seen = set()
        for g, x in zip(names, xs):
            cat = category_for(g)
            e = vol.loc[g, "sem"]
            ax.errorbar(x, x, xerr=e, yerr=e, fmt=MARKERS.get(cat, "o"),
                        color=WT_COLOR if cat == "Wild type" else COLORS[cat],
                        ecolor="0.7", elinewidth=0.9, capsize=2,
                        ms=10 if cat == "Wild type" else 7,
                        zorder=4 if cat == "Wild type" else 3,
                        label=cat if cat not in seen else None)
            seen.add(cat)
        _repel_labels(ax, fig, list(xs), list(ys), names,
                      fontsize=POINT_LABEL_FONTSIZE, n_iter=1200)
        ax.margins(0.13)
        ax.set_xlabel(r"$V$(NNIBP)  (Å$^3$)", fontsize=AXIS_LABEL_FONTSIZE)
        ax.set_ylabel(r"$V$(NNIBP)  (Å$^3$)", fontsize=AXIS_LABEL_FONTSIZE)
        ax.tick_params(axis="both", labelsize=TICK_LABEL_FONTSIZE)
        ax.annotate("measured directly in both states; no surrogate required",
                    xy=(0.5, 1.015), xycoords="axes fraction", va="bottom",
                    ha="center", fontsize=STAT_FONTSIZE)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        stats_rows.append(dict(table3_column="V(NNIBP)", surrogate="DIRECT",
                               n=len(xs), r2=1.0, p_value=0.0,
                               r2_min_leave_one_out=1.0, most_influential="",
                               spearman_rho=1.0, spearman_p=0.0))

    for ax in axes.ravel()[n:]:
        ax.set_visible(False)

    handles, labels = axes.ravel()[0].get_legend_handles_labels()
    for a in axes.ravel()[:n]:
        h2, l2 = a.get_legend_handles_labels()
        for h, l in zip(h2, l2):
            if l not in labels:
                handles.append(h); labels.append(l)
    order = [l for l in ("Wild type", "Susceptible", "Resistant", "Uncertain") if l in labels]
    fig.legend([handles[labels.index(l)] for l in order], order,
               loc="lower center", ncol=len(order), fontsize=LEGEND_FONTSIZE,
               frameon=False, bbox_to_anchor=(0.5, 0.0))

    fig.tight_layout(rect=(0, 0.06, 1, 0.985))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=args.dpi)
    fig.savefig(args.output.with_suffix(".pdf"))
    dest = args.output.parent / "tables/surrogate_validation_stats.csv"
    pd.DataFrame(stats_rows).to_csv(dest, index=False)
    logging.info("wrote %s and %s", args.output, dest)
    print(pd.DataFrame(stats_rows).to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
