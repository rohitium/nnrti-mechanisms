#!/usr/bin/env python3
"""Per-frame distributions of the eight RT-DOR contact metrics of Figure 3.

Figure 3 is a heatmap of each genotype's shift from wild type in the bound
state. It is compact and it hides shape: a contact that is bimodal, one whose
spread doubles without its mean moving, and one that is simply stable all read
as the same cell. These are the same eight quantities with the frames put
back.

Unlike the Figure 5 metrics there is no unliganded counterpart here -- every
one of these is a distance or a count between RT and doravirine, and it does
not exist without the ligand. So each genotype carries one violin, not two,
and the only reference line is wild type.

Seven of the eight come per frame from mechanism_coordinates.csv. The eighth,
the pyridinone moiety contact count, comes from the per-frame companion that
compute_dor_moiety_contacts now writes beside its per-replicate table.

Genotype order, colour and the wild-type reference follow
plot_rt_metric_violins, so the two sets of panels can be read side by side.

Usage
-----
    PYTHONPATH=src python -m nnrti.cli.plot_dor_contact_violins
"""
from __future__ import annotations

import argparse
import logging
import textwrap
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from nnrti.cli.plot_rt_metric_violins import (
    CATEGORY_ORDER,
    EXCLUDE,
    FS_LEGEND,
    FS_TICK,
    FS_TITLE,
    LW_GRID,
    LW_VIOLIN,
    LW_WT,
    MS_PUB,
    MS_REP,
    PALETTE,
    _display,
    category,
)

#: Panel label -> per-frame column in mechanism_coordinates.csv. The order is
#: Figure 3's own column order.
FRAME_METRICS = {
    "Lys103 C=O to DOR triazolinone N (Å)": "res103_bb_to_triazN",
    "Lys103 side-chain polar atoms to DOR (Å)": "res103_sc_to_dor_polar",
    "RT heavy-atom contacts (n), chlorocyanophenyl ring": "chl_ring_burial",
    "Tyr188 / chlorocyanophenyl interplanar angle (°)": "y188_interplanar_deg",
    "Residue 179 to DOR, minimum distance (Å)": "res179_to_dor_mindist",
    "Ser105 to DOR, minimum distance (Å)": "dor_to_S105_mindist",
    "Residue 106 to DOR, minimum distance (Å)": "dor_to_res106_mindist",
}

#: The one that lives in the moiety table instead.
MOIETY_METRIC = ("RT heavy-atom contacts (n), pyridinone moiety", "pyridinone")

#: Residues 103, 106 and 179 are substituted in some genotypes, so a label
#: naming the wild-type amino acid is wrong for those columns. Figure 3 says
#: "Val106" and "Val179"; Supplementary Table 4 already says "Residue 106".
#: The latter is used here, and the two Lys103 labels are kept because the
#: backbone carbonyl and the side-chain polar atoms are named by position in
#: the manuscript text as well.
NOTE_RESIDUE_NAMING = True


def load_frames(root: Path) -> dict[str, pd.DataFrame]:
    out: dict[str, pd.DataFrame] = {}
    coords = pd.read_csv(root / "results/analysis/mechanisms/mechanism_coordinates.csv")
    coords["genotype"] = coords.mutation.map(_display)
    for label, col in FRAME_METRICS.items():
        if col not in coords.columns:
            raise KeyError(f"{col} missing from mechanism_coordinates.csv")
        out[label] = pd.DataFrame({"genotype": coords.genotype,
                                   "replicate": coords.replicate,
                                   "value": coords[col]}).dropna()
    label, col = MOIETY_METRIC
    moiety = root / "results/analysis/mechanisms/dor_moiety_contacts_per_frame.csv"
    if moiety.exists():
        m = pd.read_csv(moiety)
        m["genotype"] = m.mutation.map(_display)
        out[label] = pd.DataFrame({"genotype": m.genotype, "replicate": m.replicate,
                                   "value": m[col]}).dropna()
    else:
        logging.warning("%s not found; the pyridinone panel is skipped", moiety)
    return out



def draw_contact_panel(ax, label, d, genotypes, audit=None, letter=None,
                       show_xticks=True, fs=1.0, wrap=0):
    """One contact metric on one axis.

    ``fs`` scales every font, so the standalone PNGs can carry type a reader
    sees at a glance while the eight-panel sheet, which has to fit a page,
    runs smaller.
    """
    fs_title, fs_tick = FS_TITLE * fs, FS_TICK * fs
    drawn: list[float] = []
    for i, g in enumerate(genotypes):
        v = d.loc[d.genotype == g, "value"].to_numpy(float)
        if v.size == 0:
            continue
        colour = PALETTE[category(g)][0]
        if v.size >= 10 and np.ptp(v) > 0:
            parts = ax.violinplot([v], positions=[i], widths=0.48,
                                  showextrema=False)
            for b in parts["bodies"]:
                b.set_facecolor(colour)
                b.set_edgecolor(colour)
                b.set_linewidth(LW_VIOLIN)
                b.set_alpha(0.62)
                vy = b.get_paths()[0].vertices[:, 1]
                drawn += [float(vy.min()), float(vy.max())]
        means = (d.loc[d.genotype == g].groupby("replicate").value.mean()
                 .to_numpy(float))
        ax.scatter(np.full(means.size, i), means, s=MS_PUB, facecolors="none",
                   edgecolors=colour, linewidths=1.4, zorder=5)
        ax.scatter([i], [v.mean()], s=MS_REP, c=colour, edgecolors="0.15",
                   linewidths=0.6, zorder=6)
        drawn += [float(means.min()), float(means.max()), float(v.mean())]
        if audit is not None:
            audit.append(dict(metric=label, genotype=g, n=int(v.size),
                              n_replicates=int(means.size),
                              mean=round(float(v.mean()), 3),
                              sd=round(float(v.std()), 3),
                              vmin=round(float(v.min()), 3),
                              vmax=round(float(v.max()), 3)))
    w = d.loc[d.genotype == "WT", "value"]
    if len(w):
        ax.axhline(float(w.mean()), color=PALETTE["Wild type"][0], lw=LW_WT,
                   ls="--", alpha=0.6)
    title = label if letter is None else f"{letter}.  {label}"
    # These names are long -- "RT heavy-atom contacts (n), chlorocyanophenyl
    # ring" is 50 characters -- and a left-aligned title runs off its own
    # axes into the next column. On the sheet they wrap instead.
    if wrap:
        title = textwrap.fill(title, wrap)
    ax.set_title(title, loc="left", fontsize=fs_title, pad=fs_title * 0.45)
    ax.set_xlim(-0.65, len(genotypes) - 0.35)
    # A hairline between genotypes, as on the RT-metric sheet.
    for i in range(len(genotypes) - 1):
        ax.axvline(i + 0.5, color="0.87", lw=0.7, zorder=0)
    ax.set_xticks(range(len(genotypes)))
    if show_xticks:
        ax.set_xticklabels(genotypes, rotation=55, ha="right",
                           rotation_mode="anchor", fontsize=fs_tick)
    else:
        ax.set_xticklabels([])
    lo, hi = min(drawn), max(drawn)
    pad = 0.03 * (hi - lo if hi > lo else max(abs(hi), 1.0))
    ax.set_ylim(lo - pad, hi + pad)
    ax.grid(axis="y", alpha=0.18, lw=LW_GRID)
    ax.set_axisbelow(True)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.tick_params(axis="y", labelsize=fs_tick)
    ax.tick_params(axis="x", length=fs_tick * 0.25)


def main(argv: list[str] | None = None) -> int:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    from matplotlib.legend_handler import HandlerTuple
    from matplotlib.patches import Patch

    root = Path(".")
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--manifest", type=Path, default=Path("manifests/md_manifest.csv"))
    p.add_argument("--output-dir", type=Path,
                   default=Path("results/analysis/supplementary_figures/dor_contact_violins"))
    p.add_argument("--dpi", type=int, default=200)
    p.add_argument("--no-combined", action="store_true",
                   help="Standalone panels only; skip the document sheet.")
    args = p.parse_args(argv)

    data = load_frames(root)
    man = pd.read_csv(args.manifest)[["mutation", "fold_reduction"]].drop_duplicates()
    fold = {_display(m): v for m, v in zip(man.mutation, man.fold_reduction)}
    genotypes = sorted(
        {g for d in data.values() for g in d.genotype.unique()} - EXCLUDE,
        key=lambda g: (CATEGORY_ORDER.index(category(g)) if g != "WT" else -1,
                       np.nan_to_num(fold.get(g, np.nan), nan=np.inf), g))
    logging.info("%d genotypes, %d panels", len(genotypes), len(data))

    cat_order = ["Wild type", "Susceptible", "Uncertain", "Resistant"]
    handles = [Patch(facecolor=PALETTE[c][0], alpha=0.78) for c in cat_order]
    labels = ["WT" if c == "Wild type" else c for c in cat_order]
    handles += [
        plt.Line2D([], [], marker="o", ls="", color="0.35", markersize=FS_TICK * 0.42),
        plt.Line2D([], [], marker="o", ls="", markerfacecolor="none",
                   markeredgecolor="0.35", markersize=FS_TICK * 0.5),
        plt.Line2D([], [], ls="--", color="0.35", lw=LW_WT),
    ]
    labels += ["distribution mean", "per-replicate mean", "WT"]

    args.output_dir.mkdir(parents=True, exist_ok=True)
    audit = []
    for k, (label, d) in enumerate(data.items()):
        fig, ax = plt.subplots(figsize=(13.0, 8.6))
        draw_contact_panel(ax, label, d, genotypes, audit=audit,
                           letter=chr(65 + k))
        fig.subplots_adjust(left=0.105, right=0.995, top=0.775, bottom=0.295)
        fig.legend(handles=handles, labels=labels, fontsize=FS_LEGEND * 0.72,
                   loc="lower center", ncol=3, frameon=False,
                   bbox_to_anchor=(0.55, 0.945), handlelength=2.2,
                   handler_map={tuple: HandlerTuple(ndivide=None, pad=0.0)})
        stem = FRAME_METRICS.get(label, MOIETY_METRIC[1])
        path = args.output_dir / f"{k + 1:02d}_{stem}.png"
        fig.savefig(path, dpi=args.dpi, bbox_inches="tight", pad_inches=0.25)
        plt.close(fig)
        logging.info("wrote %s", path)

    # ---- the sheet for the document ----------------------------------------
    # Same construction as the RT-metric sheet: smaller type, no panel
    # letters, one legend in a band of its own above the first row, and
    # genotype labels only under the bottom row of each column.
    if not args.no_combined:
        SHEET_FS = 0.74
        panels = list(data.items())
        nrow, ncol = 4, 2
        # An x axis under every panel, not only the bottom row -- see the
        # note in plot_rt_metric_violins.
        fig, axes = plt.subplots(nrow, ncol, figsize=(21.0, 24.5))
        for k, (label, d) in enumerate(panels):
            ax = axes[k % nrow, k // nrow]
            draw_contact_panel(ax, label, d, genotypes, audit=None,
                               letter=None, show_xticks=True,
                               fs=SHEET_FS, wrap=38)
        fig.legend(handles=handles, labels=labels, fontsize=FS_LEGEND * 0.70,
                   loc="upper center", ncol=4, frameon=False,
                   bbox_to_anchor=(0.5, 1.0), handlelength=2.6,
                   columnspacing=1.6,
                   handler_map={tuple: HandlerTuple(ndivide=None, pad=0.0)})
        fig.subplots_adjust(left=0.105, right=0.985, top=0.92, bottom=0.085,
                            hspace=0.95, wspace=0.22)
        stem = (args.output_dir.parent
                / "supplementary_figure_dor_contact_violins")
        fig.savefig(f"{stem}.png", dpi=args.dpi)
        fig.savefig(f"{stem}.pdf")
        plt.close(fig)
        logging.info("wrote %s.png / .pdf", stem)

    a = pd.DataFrame(audit)
    a.to_csv(args.output_dir / "dor_contact_violins_values.csv", index=False)
    print(a.groupby("metric").n.describe()[["count", "min", "max"]].to_string())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
