#!/usr/bin/env python3
"""Bound and unliganded distributions of all ten Figure 5 metrics, all genotypes.

Reviewer 1 asks (comment 2) for longer apo simulations of the loop-region
genotypes compared against WT, and (comment 3) whether mutation-induced
remodeling happens in both states and is therefore invisible to ddG_bind.
Figure 5 answers that as a heatmap of (holo - apo)_mut - (holo - apo)_WT: one
number and one error bar per cell. A heatmap cannot show whether a cell is a
shifted distribution, a widened one or a bimodal one, which is exactly the
distinction the reviewer's "slower conformational remodeling" turns on.

This is that figure with the distributions put back. Ten panels, one per
Figure 5 metric, nineteen genotypes on the x-axis of each, bound and
unliganded side by side.

The unit of observation differs by metric and the panel says which, because
pretending otherwise is how the first version of this figure drew violins over
three replicate means:

  volume, chi1    per frame. These exist frame by frame.
  RMSF, strand    per 10 ns block, ten per replicate. An RMSF is a fluctuation
                  about a mean structure, so it has no per-frame value and only
                  exists over a window.

A block RMSF is smaller than the 100 ns RMSF of the same replicate, because a
10 ns window samples less of the motion. The two are different quantities, so
the published 100 ns per-replicate values are overlaid as open markers rather
than left for the reader to infer from the violin. Every violin's own
per-replicate means are overlaid as filled markers.

F227C is excluded: it is the alchemical intermediate, not a clinical genotype,
and Figure 5 leaves it out too.

Genotype order and the fold-change annotations follow plot_structural_heatmap,
so the x-axis here reads in the same order as Figure 5's rows.

Usage
-----
    PYTHONPATH=src python -m nnrti.cli.plot_rt_metric_violins
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

from nnrti.analysis.panel import RESISTANT, SUSCEPTIBLE, UNCERTAIN
from nnrti.cli.plot_structural_heatmap import CATEGORY_ORDER, RT_METRICS

#: The alchemical intermediate, not a clinical genotype.
EXCLUDE = {"F227C"}

#: Phenotype category carries the hue, as in Figure 2D, and the state carries
#: the shade: darker is bound, lighter is unliganded. Colouring by state alone
#: wasted the one channel a reader uses to group genotypes, and the panel is
#: ordered by category already.
#:
#: Susceptible green and Resistant red are Figure 2D's own values. Figure 2D
#: has no wild-type category and uses grey for Uncertain; grey is given to
#: wild type here, as asked, so Uncertain moves to purple -- the colour
#: plot_pocket_distributions already used for it.
PALETTE = {
    "Wild type":   ("#000000", "#9E9E9E"),
    "Susceptible": ("#2E7D32", "#A5D6A7"),
    "Resistant":   ("#C62828", "#EF9A9A"),
    "Uncertain":   ("#6A1B9A", "#CE93D8"),
}

#: The pair offset is deliberately tight. Each genotype owns one x unit, so
#: every 0.01 taken off the offset or the violin width is white space that
#: moves from inside a genotype to between genotypes -- which is the gap a
#: reader needs at the size this sheet embeds.
STATES = (
    ("holo", "Bound (RT+DOR)", 0, -0.17),
    ("apo", "Unliganded (RT only)", 1, 0.17),
)

#: Panel label -> (source key, column). The source keys are resolved in
#: load_sources below.
SOURCES = {
    "NNIBP volume (Å³)": ("volume", "value"),
    "Tyr181 side-chain χ1 (°)": ("chi1_Tyr181", "value"),
    "Tyr188 side-chain χ1 (°)": ("chi1_Tyr188", "value"),
}

#: Crystallographic chi1 for the two pocket tyrosines, from the five
#: wild-type RT-NNRTI structures compared in Supplementary Section 1
#: (compare_nnrti_crystal_poses). Drawn on the chi1 panels as an experimental
#: anchor, because the rotamer a simulation samples is only interesting
#: against the rotamer the crystals actually show.
#:
#: The split is the point: doravirine holds Tyr181 gauche-minus in both of its
#: structures, while rilpivirine, efavirenz and nevirapine all hold it trans,
#: about 100 degrees away. A bound-state Tyr181 trans excursion in the MD is
#: therefore the side chain moving into the rotamer the OTHER NNRTIs enforce
#: and doravirine does not.
#: Tyr181 only. The five structures put Tyr188 between -46 and -69 degrees,
#: a 23 degree span that lands inside the bound violins and crowds the panel
#: without separating anything -- the informative Tyr188 result is that the
#: trans rotamer is absent from every structure AND every bound frame, which
#: is a sentence in the caption, not five more lines on the plot. Tyr181 is
#: the opposite: doravirine gauche-minus against trans for the other three
#: NNRTIs, about 100 degrees apart, and the MD visits both.
CRYSTAL_CHI1 = {
    "Tyr181": {"DOR (4NCG, 7Z2G)": (-96.0, -86.4),
               "RPV / EFV / NVP": (170.8, 168.9, 173.8)},
}
CRYSTAL_STYLE = {"DOR (4NCG, 7Z2G)": ("#111111", "-"),
                 "RPV / EFV / NVP": ("#555555", ":")}

UNIT_FRAME = "per frame"
UNIT_BLOCK = "per 10 ns block"

#: Panel titles. The quantity and its unit come first and the region after
#: the comma, so a column of panels reads down consistently. These are display
#: strings only -- RT_METRICS still carries the Figure 5 wording, so renaming a
#: panel here cannot move a figure cell or a table column.
#:
#: These are RMSF, the root-mean-square FLUCTUATION of each Ca about its own
#: mean position over the window. Not RMSD: no reference structure is involved
#: and nothing is superposed onto the crystal. The two are different
#: quantities and a reviewer would read RMSD as the latter.
#:
#: Residue ranges are deliberately absent; they belong in Methods.
DISPLAY_LABEL = {
    "NNIBP volume (\u00c5\u00b3)": "NNIBP volume (\u00c5\u00b3)",
    "Tyr181 side-chain \u03c71 (\u00b0)": "Tyr181 side-chain \u03c71 (\u00b0)",
    "Tyr188 side-chain \u03c71 (\u00b0)": "Tyr188 side-chain \u03c71 (\u00b0)",
    "C\u03b1 RMSF, NNIBP lining (\u00c5)": "C\u03b1 RMSF (\u00c5), NNIBP lining",
    "C\u03b1 RMSF, primer grip 225\u2013235 (\u00c5)": "C\u03b1 RMSF (\u00c5), primer grip",
    "C\u03b1 RMSF, YMDD motif 183\u2013186 (\u00c5)": "C\u03b1 RMSF (\u00c5), YMDD motif",
    "C\u03b1 RMSF, p66 fingers (\u00c5)": "C\u03b1 RMSF (\u00c5), p66 fingers",
    "C\u03b1 RMSF, p66 thumb (\u00c5)": "C\u03b1 RMSF (\u00c5), p66 thumb",
    "C\u03b1 RMSF, p66 connection (\u00c5)": "C\u03b1 RMSF (\u00c5), p66 connection",
    "\u03b29\u2013\u03b210 strand content (fraction)": "\u03b29\u2013\u03b210 strand content (fraction)",
}

#: Point sizes. The figure is 13 inches wide and is embedded at 6.5, so every
#: size here is halved in the document: 26 pt of tick label arrives as 13 pt,
#: about four times what the first version of this figure delivered.
FS_TITLE = 40.0
FS_SPEC = 24.0
FS_TICK = 26.0
FS_LEGEND = 27.0

LW_VIOLIN = 1.0
LW_FLAT = 3.0
LW_WT = 2.2
LW_GRID = 1.0
MS_REP = 90.0
MS_PUB = 210.0


def _display(label: str) -> str:
    return "WT" if str(label).lower() == "wt" else str(label).upper().replace("_", "+")


def category(g: str) -> str:
    if g == "WT":
        return "Wild type"
    for grp, name in ((SUSCEPTIBLE, "Susceptible"), (RESISTANT, "Resistant"),
                      (UNCERTAIN, "Uncertain")):
        if g in grp:
            return name
    return "Uncertain"


def circular_mean(deg) -> float:
    r = np.deg2rad(np.asarray(deg, dtype=float))
    return float(np.rad2deg(np.arctan2(np.sin(r).mean(), np.cos(r).mean())))


def quiet_branch_cut(values: np.ndarray, smooth_deg: int = 15) -> float:
    """The angle a chi1 axis should wrap at: the least populated one.

    A violin is a kernel density on a straight line, and chi1 lives on a
    circle. On the signed (-180, +180] convention the wrap falls at the trans
    rotamer, which is populated, so a trans side chain is drawn as two half
    violins pinned to opposite edges of the panel with the kernel smeared
    through the middle -- the one place the side chain never is. Tyr188 in the
    unliganded state is the worst case and it is not a minor cosmetic point:
    the picture asserts a bimodality that the data do not contain.

    Rotating the axis to cut at the emptiest angle fixes that without changing
    a single number. The tick labels are still written on the signed
    convention of Supplementary Table 4, so the figure and the table remain
    readable against each other.
    """
    h = np.histogram(((np.asarray(values, float) + 180.0) % 360.0) - 180.0,
                     bins=360, range=(-180.0, 180.0))[0].astype(float)
    k = np.ones(smooth_deg) / smooth_deg
    dens = np.convolve(np.r_[h, h, h], k, mode="same")[360:720]
    return float(np.argmin(dens) - 180.0)


def wrap_from(values, cut: float) -> np.ndarray:
    """Angles placed on [cut, cut + 360)."""
    return cut + (np.asarray(values, float) - cut) % 360.0


def signed(angle: float) -> float:
    return (angle + 180.0) % 360.0 - 180.0


def load_sources(root: Path, holo_volume: Path | None = None) -> dict[str, pd.DataFrame]:
    """Long-form (genotype, state, replicate, value) frames, one per metric."""
    out: dict[str, pd.DataFrame] = {}

    # ---- NNIBP volume, per frame -------------------------------------------
    apo = pd.read_csv(root / "results/analysis/apo_ensembles/tables/apo_frame_metrics.csv")
    apo = apo[apo.window == "matched"]
    if "pocket_intact" in apo.columns:
        # The same filter build_structural_metrics_table applies before it takes
        # the mean that Supplementary Table 4 prints. A disintegrated pocket has
        # an enormous convex hull; leaving those frames in put a 1871 A^3 mean in
        # circulation for apo G190E.
        apo = apo[apo.pocket_intact]
    apo = pd.DataFrame({"genotype": apo.safe_label.map(_display), "state": "apo",
                        "replicate": apo.replicate, "value": apo.pocket_volume})
    holo_path = holo_volume or (
        root / "results/analysis/pocket_volume/tables/holo_volume_per_frame.csv")
    holo = pd.read_csv(holo_path)
    holo = pd.DataFrame({"genotype": holo.mutation.map(_display), "state": "holo",
                         "replicate": holo.replicate, "value": holo.pocket_volume})

    # Three unliganded frames are dropped here, and they are the only ones in
    # the panel: G190E replicate 2, frames 39-41 at 20.0, 20.5 and 21.0 ns.
    #
    # They are the onset of the unfolding the caption already documents from
    # 21.5 ns, caught three frames too early by the pocket-Ca radius-of-
    # gyration filter -- their Rg is 8.7-9.2 A against that filter's 10 A
    # ceiling, rising but not yet over it. They are contiguous and terminal:
    # the surviving window for that replicate is frames 0-41, and these are
    # its last three.
    #
    # The test applied is physical rather than statistical. An unliganded
    # pocket cannot be larger than the largest pocket in the panel that has
    # doravirine inside it; past that the cavity has opened to bulk solvent
    # and the grid proxy is measuring solvent, not a pocket. These three
    # frames (578, 597 and 942 A^3) are the only unliganded frames in all
    # twenty genotypes above that bound, and 942 A^3 is 22% of the whole
    # 10 A scoring sphere.
    #
    # This is not a percentile clip and it removes no tail: 3 frames of
    # 14,182. Keeping them moved the G190E unliganded mean from 46.1 to
    # 102.6 A^3 on 42 surviving frames, understating the genotype's
    # induced-fit volume change as +72 +/- 26 A^3 where it is +87 +/- 16.
    ceiling = float(holo.value.max())
    n_before = len(apo)
    apo = apo[apo.value <= ceiling]
    if len(apo) != n_before:
        logging.info("volume: dropped %d unliganded frame(s) above the largest "
                     "bound pocket (%.0f A^3)", n_before - len(apo), ceiling)
    out["volume"] = pd.concat([holo, apo], ignore_index=True)

    # ---- chi1, per frame ----------------------------------------------------
    chi = pd.read_csv(root / "results/analysis/chi1_apo_holo/tables/chi1_per_frame.csv")
    chi["genotype"] = chi.mutation.map(_display)
    for residue in ("Tyr181", "Tyr188"):
        d = chi[chi.residue == residue]
        out[f"chi1_{residue}"] = pd.DataFrame({
            "genotype": d.genotype, "state": d.state,
            "replicate": d.replicate, "value": d.chi1_deg})

    # ---- RMSF and strand content, per 10 ns block ---------------------------
    blocks = pd.read_csv(root / "results/analysis/block_metrics/tables/block_metrics.csv")
    for label, col in RT_METRICS.items():
        if label in SOURCES:
            continue
        if col not in blocks.columns:
            raise KeyError(f"{col} missing from block_metrics.csv for panel {label!r}")
        out[col] = pd.DataFrame({
            "genotype": blocks.genotype, "state": blocks.state,
            "replicate": blocks.replicate, "value": blocks[col]})
        SOURCES[label] = (col, "value")
    return out


def draw_panel(ax, label, data, published, genotypes, ticks, audit,
               letter=None, show_xticks=True, fs=1.0, inside_crystal=False):
    """One metric on one axis. Returns the unit of observation it used.

    ``fs`` scales every font in the panel. The standalone PNGs keep 1.0,
    which is type a reader can see across a room; the ten-panel sheet runs
    smaller because ten panels at that size do not fit a page.

    ``inside_crystal`` moves the experimental chi1 labels inside the axes.
    They hang off the right-hand edge by default, which is fine for a lone
    panel and runs straight into the neighbouring column on the sheet.
    """
    fs_title, fs_tick = FS_TITLE * fs, FS_TICK * fs
    key, col = SOURCES[label]
    d = data[key]
    circular = "\u03c71" in label
    unit = UNIT_FRAME if (key == "volume" or key.startswith("chi1")) else UNIT_BLOCK
    pub_col = RT_METRICS[label]
    panel = d[d.genotype.isin(genotypes)]
    cut = quiet_branch_cut(panel[col].dropna().to_numpy(float)) if circular else None
    # Limits come from the artists, not from the data. A violin is a kernel
    # density and its body reaches past the sample it was built from, so an
    # axis set to the data range slices the top and bottom off the widest
    # violins -- which is what "truncated" looks like.
    drawn: list[float] = []

    for state, _state_name, shade, offset in STATES:
        for i, g in enumerate(genotypes):
            colour = PALETTE[category(g)][shade]
            v = d.loc[(d.genotype == g) & (d.state == state), col].dropna().to_numpy(float)
            if v.size == 0:
                continue
            raw = v
            if circular:
                v = wrap_from(v, cut)
            x = i + offset
            if v.size >= 10 and np.ptp(v) > 0:
                parts = ax.violinplot([v], positions=[x], widths=0.30,
                                      showextrema=False)
                for b in parts["bodies"]:
                    b.set_facecolor(colour)
                    b.set_edgecolor(colour)
                    b.set_linewidth(LW_VIOLIN)
                    b.set_alpha(0.62)
                    vy = b.get_paths()[0].vertices[:, 1]
                    drawn += [float(vy.min()), float(vy.max())]
            else:
                ax.plot([x - 0.14, x + 0.14], [np.mean(v)] * 2, color=colour,
                        lw=LW_FLAT)
                drawn += [float(np.mean(v))]
            # Open circles: the per-replicate means of the plotted
            # distribution. Filled circle: the mean of that whole
            # distribution. Both are computed from exactly what is drawn, so
            # neither can drift away from the violin it sits on.
            reps = d.loc[(d.genotype == g) & (d.state == state)]
            means = (reps.groupby("replicate")[col]
                     .apply(circular_mean if circular else "mean").to_numpy(float))
            if circular:
                means = wrap_from(means, cut)
            ax.scatter(np.full(means.size, x), means, s=MS_PUB, facecolors="none",
                       edgecolors=colour, linewidths=1.4, zorder=5)
            overall = circular_mean(raw) if circular else float(raw.mean())
            if circular:
                overall = float(wrap_from(overall, cut))
            ax.scatter([x], [overall], s=MS_REP, c=colour, edgecolors="0.15",
                       linewidths=0.6, zorder=6)
            drawn += [float(means.min()), float(means.max()), overall]
            # A metric with no published column -- the entrance gates -- has
            # nothing to overlay, which is not an error.
            if pub_col in published.columns:
                pub = published.loc[(published.genotype == g)
                                    & (published.state == state), pub_col]
                pub = pub.dropna().to_numpy(float)
            else:
                pub = np.empty(0)
            if audit is not None:
                audit.append(dict(metric=label, genotype=g, state=state, unit=unit,
                                  n=int(raw.size), n_replicates=int(means.size),
                                  mean=round(float(circular_mean(raw) if circular
                                                   else raw.mean()), 3),
                                  sd=round(float(np.std(raw)), 3),
                                  vmin=round(float(raw.min()), 3),
                                  vmax=round(float(raw.max()), 3),
                                  published_mean=(round(float(np.mean(pub)), 3)
                                                  if pub.size else np.nan)))

    shown = DISPLAY_LABEL.get(label, label)
    title = shown if letter is None else f"{letter}.  {shown}"
    ax.set_title(title, loc="left", fontsize=fs_title, pad=fs_title * 0.45)
    # On the sheet the experimental chi1 labels have to live inside the
    # axes, so the two dihedral panels get an empty gutter on the right for
    # them to sit in rather than being laid over G190E's violins.
    gutter = 2.6 if (inside_crystal and circular) else 0.0
    ax.set_xlim(-0.65, len(genotypes) - 0.35 + gutter)
    # A hairline between genotypes. At the size this sheet embeds, the eye
    # otherwise groups the apo violin of one genotype with the holo violin
    # of the next.
    for i in range(len(genotypes) - 1):
        ax.axvline(i + 0.5, color="0.87", lw=0.7, zorder=0)
    ax.set_xticks(range(len(genotypes)))
    if show_xticks:
        ax.set_xticklabels(ticks, rotation=55, ha="right",
                           rotation_mode="anchor", fontsize=fs_tick)
    else:
        ax.set_xticklabels([])

    lo, hi = min(drawn), max(drawn)
    if circular:
        pad = max(0.03 * (hi - lo), 5.0)
        # Deliberately NOT clamped to [cut, cut + 360]: clamping cut the top
        # and bottom off every violin whose density reached the wrap.
        ax.set_ylim(lo - pad, hi + pad)
        tk = [t for t in np.arange(-720.0, 1081.0, 45.0)
              if ax.get_ylim()[0] <= t <= ax.get_ylim()[1]]
        ax.set_yticks(tk)
        ax.set_yticklabels([f"{signed(t):.0f}" for t in tk])
    else:
        pad = 0.03 * (hi - lo if hi > lo else max(abs(hi), 1.0))
        ax.set_ylim(lo - pad, hi + pad)
    for state, _n, shade, _o in STATES:
        colour = PALETTE["Wild type"][shade]
        w = d.loc[(d.genotype == "WT") & (d.state == state), col].dropna()
        if len(w):
            y = (wrap_from(circular_mean(w), cut) if circular else w.mean())
            ax.axhline(float(y), color=colour, lw=LW_WT, ls="--", alpha=0.6)
    if circular:
        residue = "Tyr181" if "Tyr181" in label else "Tyr188"
        lo_a, hi_a = ax.get_ylim()
        for name, values in CRYSTAL_CHI1.get(residue, {}).items():
            colour, style = CRYSTAL_STYLE[name]
            ys = [float(wrap_from(v, cut)) for v in values]
            ys = [y for y in ys if lo_a <= y <= hi_a]
            if not ys:
                continue
            for y in ys:
                ax.axhline(y, color=colour, ls=style, lw=LW_WT * 0.85,
                           alpha=0.85, zorder=1)
            if inside_crystal:
                ax.annotate(name, xy=(len(genotypes) - 0.15,
                                      float(np.mean(ys))),
                            xytext=(0, fs_tick * 0.25),
                            textcoords="offset points",
                            va="bottom", ha="left", fontsize=fs_tick * 0.78,
                            color=colour)
            else:
                ax.annotate(name, xy=(1.0, float(np.mean(ys))),
                            xycoords=("axes fraction", "data"),
                            xytext=(fs_tick * 0.25, 0),
                            textcoords="offset points",
                            va="center", ha="left", fontsize=fs_tick * 0.80,
                            color=colour)

    ax.grid(axis="y", alpha=0.18, lw=LW_GRID)
    ax.set_axisbelow(True)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.tick_params(axis="y", labelsize=fs_tick)
    ax.tick_params(axis="x", length=fs_tick * 0.25)
    return unit


def legend_handles(fs_tick):
    from matplotlib.patches import Patch
    cat_order = ["Wild type", "Susceptible", "Uncertain", "Resistant"]
    handles = [(Patch(facecolor=PALETTE[c][0], alpha=0.78),
                Patch(facecolor=PALETTE[c][1], alpha=0.78)) for c in cat_order]
    labels = [f"{'WT' if c == 'Wild type' else c}: bound, unliganded"
              for c in cat_order]
    handles += [
        plt.Line2D([], [], marker="o", ls="", color="0.35", markersize=fs_tick * 0.42),
        plt.Line2D([], [], marker="o", ls="", markerfacecolor="none",
                   markeredgecolor="0.35", markersize=fs_tick * 0.5),
        plt.Line2D([], [], ls="--", color="0.35", lw=LW_WT),
    ]
    labels += ["distribution mean", "per-replicate mean", "WT, same state"]
    return handles, labels


#: Predicted resistant if the FEP ddG_bind, less its standard error, clears
#: the manuscript's 0.5 kcal/mol "strong impact" threshold -- the same rule
#: classification_performance.py applies. The reviewer's three mismatch cases
#: are named against the FEP numbers, not the MM/GBSA ones: K103N+M230L
#: +0.59 +/- 0.41 at 36x, A98G+F227C -1.53 +/- 0.87 at 93x, G190E
#: -1.02 +/- 0.38 at 18x.
DDG_THRESHOLD = 0.5


def prediction_buckets(fep_csv: Path) -> dict[str, str]:
    """genotype -> 'Predicted susceptible' | 'Predicted resistant'."""
    f = pd.read_csv(fep_csv)
    out = {}
    for _, r in f.iterrows():
        g = _display(r["genotype"])
        pred = (r["ddg_bind_kcal"] - r["sem_kcal"]) > DDG_THRESHOLD
        out[g] = "Predicted resistant" if pred else "Predicted susceptible"
    return out


def draw_bucketed(label, data, published, buckets, fold, out_dir, letter, dpi):
    """One metric, the two prediction buckets side by side on a shared y-axis.

    The question this answers is not whether resistance separates from
    susceptibility -- the panel order already carries that -- but whether the
    genotypes ddG_bind gets wrong look different from the ones it gets right.
    So the split is by what the cycle PREDICTED, and the colour inside each
    half stays the observed phenotype: a red violin in the left-hand panel is
    a genotype the calculation called harmless and the clinic did not.
    """
    import matplotlib.pyplot as plt
    from matplotlib.legend_handler import HandlerTuple

    groups = ["Predicted susceptible", "Predicted resistant"]
    members = {b: sorted([g for g, v in buckets.items() if v == b],
                         key=lambda g: np.nan_to_num(fold.get(g, np.nan), nan=0.0))
               for b in groups}
    widths = [len(members[b]) for b in groups]
    fig, axes = plt.subplots(
        1, 2, figsize=(15.0, 9.0), sharey=True,
        gridspec_kw={"width_ratios": [max(w, 1) for w in widths]})
    for ax, b in zip(axes, groups):
        gl = ["WT"] + members[b]
        draw_panel(ax, label, data, published, gl, list(gl), None,
                   letter=None, show_xticks=True)
        ax.set_title(b, loc="left", fontsize=FS_TITLE * 0.62,
                     pad=FS_TITLE * 0.25)
    axes[1].tick_params(axis="y", left=False)
    shown = DISPLAY_LABEL.get(label, label)
    fig.suptitle(f"{letter}.  {shown}", x=0.015, ha="left",
                 fontsize=FS_TITLE, y=0.995)
    handles, labels = legend_handles(FS_TICK)
    # Title band on top, legend under it, axes under that: the three were
    # overlapping when the legend was anchored near the top of the canvas.
    fig.subplots_adjust(left=0.095, right=0.995, top=0.675, bottom=0.30,
                        wspace=0.04)
    fig.legend(handles=handles, labels=labels, fontsize=FS_LEGEND * 0.66,
               loc="lower center", ncol=2, frameon=False,
               bbox_to_anchor=(0.56, 0.745), handlelength=3.0,
               handler_map={tuple: HandlerTuple(ndivide=None, pad=0.0)})
    path = out_dir / f"{letter}_{SOURCES[label][0]}_by_prediction.png"
    fig.savefig(path, dpi=dpi, bbox_inches="tight", pad_inches=0.25)
    plt.close(fig)
    return path


def main(argv: list[str] | None = None) -> int:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    from matplotlib.legend_handler import HandlerTuple

    root = Path(".")
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--manifest", type=Path, default=Path("manifests/md_manifest.csv"))
    p.add_argument("--published", type=Path,
                   default=Path("results/analysis/structural_metrics/tables/"
                                "structural_metrics_per_replicate.csv"))
    p.add_argument("--holo-volume", type=Path, default=None)
    p.add_argument("--output-dir", type=Path,
                   default=Path("results/analysis/supplementary_figures"))
    p.add_argument("--no-combined", action="store_true",
                   help="Skip the stacked all-panel sheet.")
    p.add_argument("--gates", action="store_true",
                   help="Also draw the NNIBP entrance-gate distances.")
    p.add_argument("--gate-csv", type=Path,
                   default=Path("results/analysis/nnibp_gates/tables/"
                                "nnibp_gate_distances_per_frame.csv"))
    p.add_argument("--by-prediction", action="store_true",
                   help="Also split each metric by predicted ddG_bind bucket.")
    p.add_argument("--fep", type=Path,
                   default=Path("results/analysis/fep_pmx/panel_ddg.csv"))
    p.add_argument("--dpi", type=int, default=200)
    args = p.parse_args(argv)

    data = load_sources(root, args.holo_volume)
    published = pd.read_csv(args.published).rename(columns={"structure_type": "state"})

    man = pd.read_csv(args.manifest)[["mutation", "fold_reduction"]].drop_duplicates()
    fold = {_display(m): v for m, v in zip(man.mutation, man.fold_reduction)}

    genotypes = sorted(
        {g for d in data.values() for g in d.genotype.unique()} - EXCLUDE,
        key=lambda g: (CATEGORY_ORDER.index(category(g)) if g != "WT" else -1,
                       np.nan_to_num(fold.get(g, np.nan), nan=np.inf), g))
    # Genotype alone. The fold reductions are in Table 1 and repeating them on
    # every one of ten panels crowded the axis without adding anything the
    # category colour does not already carry.
    ticks = list(genotypes)
    metrics = list(RT_METRICS)
    logging.info("%d genotypes, %d panels", len(genotypes), len(metrics))

    out = args.output_dir / "rt_metric_violins"
    out.mkdir(parents=True, exist_ok=True)
    audit: list[dict] = []
    handles, labels = legend_handles(FS_TICK)

    # ---- one standalone PNG per metric, each with its own genotype axis ----
    for k, label in enumerate(metrics):
        fig, ax = plt.subplots(figsize=(13.0, 8.6))
        draw_panel(ax, label, data, published, genotypes, ticks, audit,
                   letter=chr(65 + k), show_xticks=True)
        # The legend sits directly above the axes, not at the top of a blank
        # figure: anchoring it to 1.0 while the axes stop at 0.7 leaves a band
        # of white between the two that bbox_inches="tight" cannot remove,
        # because the legend is what defines the top of the tight box.
        fig.subplots_adjust(left=0.105, right=0.80, top=0.775, bottom=0.295)
        fig.legend(handles=handles, labels=labels, fontsize=FS_LEGEND * 0.72,
                   loc="lower center", ncol=2, frameon=False,
                   bbox_to_anchor=(0.55, 0.955), handlelength=3.0,
                   handler_map={tuple: HandlerTuple(ndivide=None, pad=0.0)})
        name = SOURCES[label][0].replace("chi1_", "chi1_").replace("/", "_")
        path = out / f"{k + 1:02d}_{name}.png"
        fig.savefig(path, dpi=args.dpi, bbox_inches="tight", pad_inches=0.25)
        plt.close(fig)
        logging.info("wrote %s", path)

    # ---- the stacked sheet, for the document --------------------------------
    if not args.no_combined:
        # Two columns for the document copy. One column at the standalone
        # panels' font size is 54 inches tall, which embeds at 27 inches and
        # cannot fit a page -- Word does not break an image across pages, so
        # it would simply overflow. The standalone PNGs keep the large type;
        # this sheet trades some of it for a figure that fits.
        nrow, ncol = 5, 2
        # SHEET_FS shrinks every font in the panels. At 1.0 the titles of the
        # right-hand column ran off the canvas, the rotated genotype labels
        # of the left-hand column ran off the left of it, and the legend band
        # sat on top of the first row of titles. Panel letters are gone too:
        # ten of them label nothing a caption cannot say in order.
        SHEET_FS = 0.74
        # Every panel carries its own genotype axis. Labelling only the
        # bottom row saved vertical space but made the upper eight panels
        # unreadable: a reader had to count violins down four panels to
        # learn which genotype a feature belonged to.
        fig, axes = plt.subplots(nrow, ncol, figsize=(21.0, 27.8))
        for k, label in enumerate(metrics):
            ax = axes[k % nrow, k // nrow]
            draw_panel(ax, label, data, published, genotypes, ticks, None,
                       letter=None, show_xticks=True,
                       fs=SHEET_FS, inside_crystal=True)
        fig.legend(handles=handles, labels=labels, fontsize=FS_LEGEND * 0.70,
                   loc="upper center", ncol=4, frameon=False,
                   bbox_to_anchor=(0.5, 1.0), handlelength=2.6,
                   columnspacing=1.6,
                   handler_map={tuple: HandlerTuple(ndivide=None, pad=0.0)})
        # top leaves a clear band for the two legend rows; left and bottom
        # leave room for the rotated genotype labels; right stops the second
        # column's title at the canvas edge.
        # bottom has to clear a full set of rotated genotype labels now that
        # the last row carries one too; left clears the longest of them,
        # which overhangs the first panel's y axis.
        # 21 x 27.8 in is deliberate: at that aspect the sheet embeds at the
        # full 6.5 in text width without exceeding the page height, which is
        # what sets how far apart the genotypes land on paper.
        fig.subplots_adjust(left=0.105, right=0.985, top=0.935, bottom=0.062,
                            hspace=1.00, wspace=0.22)
        stem = args.output_dir / "supplementary_figure_rt_metric_violins"
        fig.savefig(f"{stem}.png", dpi=args.dpi)
        fig.savefig(f"{stem}.pdf")
        plt.close(fig)
        logging.info("wrote %s.png / .pdf", stem)

    if args.gates and args.gate_csv.exists():
        # The gates are RT-only two-state metrics, exactly the shape
        # draw_panel already handles, so they are registered into the same
        # tables rather than given a parallel drawing path.
        gates = pd.read_csv(args.gate_csv)
        gdir = args.output_dir / "nnibp_gate_violins"
        gdir.mkdir(parents=True, exist_ok=True)
        for n, (gate, sub) in enumerate(gates.groupby("gate")):
            label = f"NNIBP entrance, {gate} minimum distance (\u00c5)"
            key = "gate_" + gate.split()[0].replace("-", "_")
            data[key] = pd.DataFrame({"genotype": sub.genotype, "state": sub.state,
                                      "replicate": sub.replicate,
                                      "value": sub.min_dist_A})
            SOURCES[label] = (key, "value")
            RT_METRICS[label] = "__gate__"   # no published column to overlay
            fig, ax = plt.subplots(figsize=(13.0, 8.6))
            draw_panel(ax, label, data, published, genotypes, ticks, None,
                       letter=None, show_xticks=True)
            fig.subplots_adjust(left=0.105, right=0.995, top=0.775, bottom=0.295)
            fig.legend(handles=handles, labels=labels, fontsize=FS_LEGEND * 0.72,
                       loc="lower center", ncol=2, frameon=False,
                       bbox_to_anchor=(0.55, 0.955), handlelength=3.0,
                       handler_map={tuple: HandlerTuple(ndivide=None, pad=0.0)})
            path = gdir / f"{n + 1:02d}_{key}.png"
            fig.savefig(path, dpi=args.dpi, bbox_inches="tight", pad_inches=0.25)
            plt.close(fig)
            logging.info("wrote %s", path)

    if args.by_prediction:
        buckets = {g: b for g, b in prediction_buckets(args.fep).items()
                   if g in genotypes and g != "WT"}
        bdir = args.output_dir / "rt_metric_violins_by_prediction"
        bdir.mkdir(parents=True, exist_ok=True)
        for k, label in enumerate(metrics):
            path = draw_bucketed(label, data, published, buckets, fold, bdir,
                                 chr(65 + k), args.dpi)
            logging.info("wrote %s", path)

    a = pd.DataFrame(audit)
    a.to_csv(args.output_dir / "supplementary_figure_rt_metric_violins_values.csv",
             index=False)
    print(a.groupby(["unit", "state"]).n.describe()[["count", "min", "max"]].to_string())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
