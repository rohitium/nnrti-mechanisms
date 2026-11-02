#!/usr/bin/env python3
"""Figure 3: structural change across the panel, in two heatmaps.

Reviewer 1 (comment 4) asked for a WT-referenced heatmap summarising the
structural metrics, and (comment 3) for whether mutation-induced changes occur
in both apo and holo RT -- because a mutation that moves both states equally is
invisible to ddG_bind by construction.

The two panels answer different questions and so are scaled differently:

  A  DOR-RT contacts, holo, as a shift from wild type. These need a ligand, so
     they exist only in the bound state and the apo comparison is undefined.

  B  RT's own conformation as delta = holo - apo, per genotype. This is the
     induced fit: how much binding DOR changes the enzyme. Wild type is kept
     as the first row, because the question is not whether a genotype has an
     induced fit but whether its induced fit differs from wild type's. A row
     that matches WT is a mutation whose effect on RT is concordant between
     the two states, and therefore cancels in ddG_bind.

Cells carry the value in the metric's own units above its standard error, with
errors added in quadrature wherever two means are subtracted. Colour is scaled
within each column and carries no meaning across columns; the number is the
quantity.

Usage
-----
    PYTHONPATH=src python -m nnrti.cli.plot_structural_heatmap
"""

from __future__ import annotations

import argparse
import textwrap
from pathlib import Path

import numpy as np
import pandas as pd

#: Panel A: contacts between RT and DOR. Undefined without a ligand.
#:
#: Names are copied verbatim from Table 3 of the manuscript, units included,
#: so a reader moving between the table and the figure meets one wording for
#: each quantity. The compact symbolic forms (d(103 C=O...N4x) and the like)
#: live only in the main-text table's own header row.
DOR_METRICS = {
    "Lys103 C=O to DOR triazolinone N (\u00c5)": "res103_CO_to_triazolinone_N_A",
    "Lys103 side-chain polar atoms to DOR (\u00c5)": "res103_sidechain_to_DOR_A",
    "RT heavy-atom contacts < 4.0 \u00c5 of chlorocyanophenyl ring":
        "chlorocyanophenyl_contacts_n",
    "Tyr188 / chlorocyanophenyl interplanar angle (\u00b0)":
        "Tyr188_chlorophenyl_angle_deg",
    "Val179 to DOR minimum distance (\u00c5)": "res179_to_DOR_A",
    "Ser105 to DOR minimum distance (\u00c5)": "Ser105_to_DOR_A",
    "Val106 to DOR minimum distance (\u00c5)": "res106_to_DOR_A",
    "Pyridinone moiety contacts < 4.0 \u00c5": "pyridinone_contacts_n",
}

#: Panel B: RT's own conformation, measurable with or without the ligand.
#: Named in the same style, since these are not in Table 3.
RT_METRICS = {
    "NNIBP volume (\u00c5\u00b3)": "nnibp_volume_A3",
    "Tyr181 side-chain \u03c71 (\u00b0)": "chi1_Tyr181_deg",
    "Tyr188 side-chain \u03c71 (\u00b0)": "chi1_Tyr188_deg",
    "C\u03b1 RMSF, NNIBP lining (\u00c5)": "rmsf_A_nnibp_lining",
    "C\u03b1 RMSF, primer grip 225\u2013235 (\u00c5)": "rmsf_A_primer_grip",
    "C\u03b1 RMSF, YMDD motif 183\u2013186 (\u00c5)": "rmsf_A_ymdd",
    "C\u03b1 RMSF, p66 fingers (\u00c5)": "rmsf_A_fingers",
    "C\u03b1 RMSF, p66 thumb (\u00c5)": "rmsf_A_thumb",
    "C\u03b1 RMSF, p66 connection (\u00c5)": "rmsf_A_connection",
    "\u03b29\u2013\u03b210 strand content (fraction)": "strand_frac_beta9_beta10",
}

#: Dihedrals: mean and difference must both be taken on the circle.
#: Dihedral metrics, keyed off RT_METRICS so a rename cannot silently turn the
#: circular handling off. It did exactly that once: the labels were rewritten
#: to the Table 3 wording while this set still held the old names, nothing
#: matched, and Y318F's Tyr188 shift printed as -272 degrees instead of +88.
CIRCULAR = {k for k, v in RT_METRICS.items() if v.startswith("chi1_")}
assert CIRCULAR, "no dihedral metrics found; check the chi1_ column prefix"


#: Points to the left of the axis for the category label. The longest tick
#: label ("V106I+F227C  105x", 17 characters at 10 pt) runs to about 102 pt, and
#: the rotated category text is ~12 pt wide about its centre, so anything under
#: ~115 pt overlaps and clips the first letter off the longest genotypes.
#: Clearance in points for the category label. It has to exceed the widest
#: genotype tick label, which scales with the tick font: "V106A+P225H  153x"
#: is 17 characters, so at 32 pt it runs to roughly 330 pt, and the rotated
#: category text adds its own height about its centre.
CATEGORY_LABEL_PAD_PT = 400.0

#: Standard errors a shift must exceed for its 95% interval to exclude zero.
#: The two panels differ because they are built from different numbers of
#: groups, and using one value for both is wrong in a way that matters:
#:
#:   A  mutant minus WT, two groups of three replicates -> Welch df ~ 4,
#:      t(0.975) = 2.776.
#:   B  (holo - apo)_mut minus (holo - apo)_WT, four groups of three ->
#:      df ~ 8, t(0.975) = 2.306.
#:
#: Applying A's 2.776 to panel B dropped G190E's +71 +/- 26 A^3 pocket
#: expansion, which is 2.73 SEM and comfortably significant on 8 degrees of
#: freedom. Twice the standard error, the obvious choice and the one used here
#: first, is not a 95% criterion at either size: at df = 4 it is p = 0.12.
#: Exact per-comparison p values are in structural_effect_sizes.csv.
SIGNIFICANCE_BY_PANEL = {"A": 2.776, "B": 2.306}

#: Category order of Table 1 in the manuscript.
CATEGORY_ORDER = ("Susceptible", "Resistant", "Uncertain")


def _style():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({
        "font.size": 34, "axes.labelsize": 34, "axes.titlesize": 44,
        "xtick.labelsize": 32, "ytick.labelsize": 32,
        "figure.dpi": 300, "savefig.bbox": "tight",
    })
    return plt


def _circ_mean_sem(deg: np.ndarray) -> tuple[float, float]:
    r = np.deg2rad(np.asarray(deg, float))
    r = r[np.isfinite(r)]
    if r.size == 0:
        return np.nan, np.nan
    s, c = np.sin(r).mean(), np.cos(r).mean()
    R = float(np.hypot(s, c))
    sd = float(np.rad2deg(np.sqrt(-2.0 * np.log(max(R, 1e-12)))))
    return float(np.rad2deg(np.arctan2(s, c))), sd / np.sqrt(r.size)


def _stats(d: pd.DataFrame, col: str, circular: bool) -> pd.DataFrame:
    if circular:
        out = d.groupby("genotype")[col].apply(
            lambda v: pd.Series(_circ_mean_sem(v.to_numpy()), index=["mean", "sem"]))
        return out.unstack()
    g = d.groupby("genotype")[col]
    return pd.DataFrame({"mean": g.mean(), "sem": g.sem()})


def _wrap(delta: pd.Series, circular: bool) -> pd.Series:
    return ((delta + 180.0) % 360.0) - 180.0 if circular else delta


def shift_vs_wt(per_rep, state, metrics):
    """Each genotype's shift from WT, in the metric's own units."""
    d = per_rep[per_rep.structure_type == state]
    val, err = {}, {}
    for label, col in metrics.items():
        if col not in d or d[col].isna().all():
            continue
        s = _stats(d, col, label in CIRCULAR)
        if "WT" not in s.index or not np.isfinite(s.loc["WT", "mean"]):
            continue
        val[label] = _wrap(s["mean"] - s.loc["WT", "mean"], label in CIRCULAR)
        err[label] = np.sqrt(s["sem"] ** 2 + float(s.loc["WT", "sem"]) ** 2)
    return pd.DataFrame(val), pd.DataFrame(err)


def induced_fit_vs_wt(per_rep, metrics):
    """How a genotype's induced fit differs from wild type's.

    The induced fit is delta = holo - apo: how far binding DOR moves a metric.
    Referencing that to wild type gives (holo - apo)_mut - (holo - apo)_WT, a
    difference of differences. Zero means the mutation shifts the bound and
    unliganded states by the same amount -- concordant remodelling, which
    cancels in ddG_bind by construction and is exactly the case Reviewer 1
    asks about. Four means are subtracted, so four SEMs add in quadrature.
    """
    val, err = {}, {}
    for label, col in metrics.items():
        if col not in per_rep or per_rep[col].isna().all():
            continue
        circ = label in CIRCULAR
        h = _stats(per_rep[per_rep.structure_type == "holo"], col, circ)
        a = _stats(per_rep[per_rep.structure_type == "apo"], col, circ)
        common = h.index.intersection(a.index)
        if "WT" not in common:
            continue
        fit = _wrap(h.loc[common, "mean"] - a.loc[common, "mean"], circ)
        fit_sem = np.sqrt(h.loc[common, "sem"] ** 2 + a.loc[common, "sem"] ** 2)
        val[label] = _wrap(fit - fit["WT"], circ)
        err[label] = np.sqrt(fit_sem ** 2 + float(fit_sem["WT"]) ** 2)
    return pd.DataFrame(val), pd.DataFrame(err)


def _decimals(largest: float) -> int:
    return 0 if largest >= 10 else (1 if largest >= 1 else 2)


def draw(ax, val, err, title, ylabels, show_y):
    scaled = val.astype(float).copy()
    for c in scaled.columns:
        m = float(np.nanmax(np.abs(scaled[c].to_numpy())))
        scaled[c] = scaled[c] / m if m else 0.0
    im = ax.imshow(scaled.to_numpy(dtype=float), cmap="RdBu_r", vmin=-1, vmax=1,
                   aspect="auto")

    dec = {c: _decimals(float(np.nanmax(np.abs(val[c].to_numpy()))))
           for c in val.columns}
    for i in range(val.shape[0]):
        for j, c in enumerate(val.columns):
            v = val.iat[i, j]
            if not np.isfinite(v):
                ax.text(j, i, "·", ha="center", va="center", color="0.6", fontsize=32)
                continue
            colour = "white" if abs(float(scaled.iat[i, j])) > 0.65 else "0.12"
            ax.text(j, i - 0.15, f"{v:+.{dec[c]}f}", ha="center", va="center",
                    fontsize=30, color=colour)
            e = err.iat[i, j]
            if np.isfinite(e):
                ax.text(j, i + 0.24, f"±{e:.{dec[c]}f}", ha="center", va="center",
                        fontsize=23, color=colour, alpha=0.9)

    ax.set_xticks(range(val.shape[1]))
    # Wrapped: the Table 3 names are long, and a single 45-degree line of them
    # would need more bottom margin than the heatmap itself.
    ax.set_xticklabels([textwrap.fill(str(c), 26) for c in val.columns],
                       rotation=45, ha="right", rotation_mode="anchor",
                       fontsize=26)
    ax.set_yticks(range(val.shape[0]))
    ax.set_yticklabels(ylabels if show_y else [], fontsize=32)
    ax.set_title(title, fontsize=44, pad=32)
    ax.set_xticks(np.arange(-0.5, val.shape[1], 1), minor=True)
    ax.set_yticks(np.arange(-0.5, val.shape[0], 1), minor=True)
    ax.grid(which="minor", color="white", lw=1.0)
    ax.tick_params(which="minor", length=0)
    return im


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--per-replicate", type=Path,
                   default=Path("results/analysis/structural_metrics/tables/"
                                "structural_metrics_per_replicate.csv"))
    p.add_argument("--table3-csv", type=Path,
                   default=Path("paper/tables/Table-3-structural.csv"),
                   help="Source of the Table 1 resistance categories.")
    p.add_argument("--manifest", type=Path, default=Path("manifests/md_manifest.csv"))
    p.add_argument("--output-dir", type=Path,
                   default=Path("results/analysis/structural_metrics"))
    args = p.parse_args(argv)

    plt = _style()
    from matplotlib.transforms import blended_transform_factory
    per_rep = pd.read_csv(args.per_replicate)
    man = pd.read_csv(args.manifest)[["mutation", "fold_reduction"]].drop_duplicates()
    fold = dict(zip(man.mutation, man.fold_reduction))
    cat = dict(pd.read_csv(args.table3_csv)[["Genotype", "Category"]].values)


    dor_v, dor_e = shift_vs_wt(per_rep, "holo", DOR_METRICS)
    rt_v, rt_e = induced_fit_vs_wt(per_rep, RT_METRICS)

    # Genotypes absent from Table 1 are not part of the published panel and are
    # not plotted. F227C is the only one; the build says so rather than
    # dropping it silently.
    uncategorised = sorted({g for g in per_rep.genotype.unique()
                            if g not in cat and g != "WT"})
    if uncategorised:
        print(f"dropped, not in Table 1: {', '.join(uncategorised)}")

    def key(g):
        rank = CATEGORY_ORDER.index(cat[g])
        return (rank, np.nan_to_num(fold.get(g, np.nan), nan=np.inf))

    mutants = sorted((g for g in dor_v.index if g != "WT" and g in cat), key=key)
    labels = [f"{g}  {fold.get(g):g}\u00d7" if np.isfinite(fold.get(g, np.nan))
              else g for g in mutants]

    def rule_and_label(ax, rows):
        cats = [cat[g] for g in rows]
        for i in range(1, len(cats)):
            if cats[i] != cats[i - 1]:
                ax.axhline(i - 0.5, color="black", lw=2.4)
        seen, start_i = cats[0], 0
        for i in range(1, len(cats) + 1):
            if i == len(cats) or cats[i] != seen:
                # Offset in POINTS, not axes fraction: a fraction of the
                # axes width is ~5 inches on a 20-inch-wide panel, which
                # stranded these labels far out in the margin.
                ax.annotate(seen, xy=(0, (start_i + i - 1) / 2.0),
                            xycoords=("axes fraction", "data"),
                            xytext=(-CATEGORY_LABEL_PAD_PT, 0),
                            textcoords="offset points", rotation=90,
                            ha="center", va="center", fontsize=34,
                            fontweight="bold", color="0.25", clip_on=False)
                if i < len(cats):
                    seen, start_i = cats[i], i

    out = args.output_dir
    (out / "plots").mkdir(parents=True, exist_ok=True)
    (out / "tables").mkdir(parents=True, exist_ok=True)

    def keep_resolved(val, err, rows, tag):
        """Drop metrics where not one genotype clears twice its own error.

        A column of numbers that are all smaller than their uncertainty says
        nothing about any genotype, and printing it invites the eye to read
        structure into noise. Which metrics these are is itself a result, so
        they are named rather than quietly removed.
        """
        v, e = val.loc[rows], err.reindex(rows)
        crit = SIGNIFICANCE_BY_PANEL[tag]
        sig = (v.abs() > crit * e).sum()
        drop = [c for c in v.columns if sig[c] == 0]
        if drop:
            print(f"{tag}: dropped, no genotype's 95% interval excludes zero "
                  f"({crit:.3g}x SEM) -- {', '.join(drop)}")
        keep = [c for c in v.columns if sig[c] > 0]
        return val[keep], err[keep]

    dor_v, dor_e = keep_resolved(dor_v, dor_e, mutants, "A")
    rt_v, rt_e = keep_resolved(rt_v, rt_e, mutants, "B")

    panels = [
        ("A", dor_v, dor_e, "RT\u2013DOR contacts (holo): shift from wild type",
         "figure3A_rt_dor_contacts", "figure3A_dor"),
        ("B", rt_v, rt_e,
         "RT conformation: induced fit (holo \u2212 apo), shift from wild type",
         "figure3B_rt_induced_fit", "figure3B_rt"),
    ]
    written = []
    for tag, val, err, title, stem, tstem in panels:
        v = val.loc[mutants]
        e = err.reindex(mutants)
        # Canvas deliberately NOT scaled with the font. Growing both together
        # leaves the type the same size relative to the cells and only makes
        # the file enormous -- doubling the fonts and the figure gave a 46 inch
        # wide image that read exactly as before.
        # Shrinking the canvas while the fonts stay fixed in points makes the
        # type larger relative to the cells, which is the goal, and keeps the
        # render under 8000 px on its long side so it can be viewed anywhere.
        fig, ax = plt.subplots(figsize=(2.05 * v.shape[1] + 5.0,
                                        0.92 * len(mutants) + 5.0))
        im = draw(ax, v, e, title, labels, True)
        rule_and_label(ax, mutants)
        # Colour is normalised per column, so the bar is a relative scale: the
        # darkest cell in a column is that metric's largest shift. One bar in
        # real units would be a lie, since the columns are in A, degrees, A^3
        # and counts.
        cb = fig.colorbar(im, ax=ax, fraction=0.030, pad=0.02)
        cb.set_ticks([-1, -0.5, 0, 0.5, 1])
        cb.set_ticklabels(["\u2212max", "", "0", "", "+max"])
        cb.set_label("shift, scaled within each column\n"
                     "(max = largest shift in that metric)", fontsize=30)
        cb.ax.tick_params(labelsize=30)
        f = out / "plots" / f"{stem}.png"
        fig.savefig(f, dpi=300, bbox_inches="tight", pad_inches=0.25)
        plt.close(fig)
        written.append(f)
        v.round(3).to_csv(out / "tables" / f"{tstem}_shift_vs_wt.csv")
        e.round(3).to_csv(out / "tables" / f"{tstem}_shift_sem.csv")
        n = int(np.isfinite(v.to_numpy()).sum())
        sig = int((v.abs() > SIGNIFICANCE_BY_PANEL[tag] * e).to_numpy().sum())
        print(f"{tag}: {v.shape[0]} genotypes x {v.shape[1]} metrics, "
              f"{sig}/{n} cells with a 95% interval excluding zero")

    for f in written:
        print(f"wrote {f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
