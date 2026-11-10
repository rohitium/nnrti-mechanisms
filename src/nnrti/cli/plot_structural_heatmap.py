#!/usr/bin/env python3
"""Build the current structural heatmaps from replicate means.

Main Figure 3 shows bound RT–DOR contacts, main Figure 5 compares RT metrics
with WT separately in each state, and Supplementary Figure 7 shows binding-
associated geometric changes. These geometric differences do not measure
cancellation of binding free energies. Confidence intervals and minimum
interpretable differences jointly determine highlighting.
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
#:   A  mutant minus WT, two groups of three replicates -> pooled df = 4,
#:      t(0.975) = 2.776.
#:   B  (holo - apo)_mut minus (holo - apo)_WT, four groups of three ->
#:      df = 8, t(0.975) = 2.306.
#:
#: Applying A's 2.776 to panel B dropped G190E's +71 +/- 26 A^3 pocket
#: expansion, which is 2.73 SEM and comfortably significant on 8 degrees of
#: freedom. Twice the standard error, the obvious choice and the one used here
#: first, is not a 95% criterion at either size: at df = 4 it is p = 0.12.
#: Exact per-comparison p values are in structural_effect_sizes.csv.
#: Panels C and D are mutant-minus-WT within one state: two groups of three
#: replicates, pooled df = 4, the same criterion as panel A.
SIGNIFICANCE_BY_PANEL = {"A": 2.776, "B": 2.306, "C": 2.776, "D": 2.776}

#: Smallest difference worth reading on each observable -- the "minimum
#: interpretable difference". A t test asks only whether a shift is
#: distinguishable from zero; it never asks whether the shift is large enough
#: to mean anything. With three replicates those two questions come apart: a
#: genotype whose replicates happen to agree closely can clear the 95%
#: criterion on a shift that is chemically nothing. K103N+P225H did exactly
#: that on the Lys103 hydrogen bond, clearing it on +0.13 A because its three
#: replicates spanned 0.08 A, while V106I's +0.65 A -- a real loss of the bond
#: in two replicates of three -- did not. A cell is now filled only when its
#: 95% interval excludes zero AND its magnitude clears the floor below.
#:
#: The floors are set from what the observable measures, not from which cells
#: they remove, and each applies to every panel:
#:
#:   0.30 A   interatomic distances. A hydrogen bond at 3.10 A and one at
#:            2.97 A are the same hydrogen bond.
#:   1.0      contact counts, which are integer-valued per frame.
#:   5 deg    the interplanar angle and the chi1 dihedrals, well inside one
#:            rotamer well.
#:   10 A^3   NNIBP volume, ~4% of the holo wild-type pocket.
#:   0.05 A   Ca RMSF, which is a fluctuation amplitude of a few tenths of an
#:            Angstrom, so the distance floor would be the whole quantity.
#:   0.05     secondary-structure fractions, i.e. five percentage points.
#:
#: A check that these are not reverse-engineered: applied to the Val106-DOR
#: minimum distance the 0.30 A floor independently removes both of that
#: column's resolved cells (+0.27 and +0.28 A), which is the judgement that
#: had already taken the column out of the drawn panel by hand.
MIN_EFFECT_BY_UNIT = {
    "\u00c5": 0.30,
    "contacts": 1.0,
    "degrees": 5.0,
    "\u00c5\u00b3": 10.0,
    "fraction": 0.05,
}

#: Per-metric overrides where the unit alone sets the wrong scale.
MIN_EFFECT = {
    "C\u03b1 RMSF, NNIBP lining (\u00c5)": 0.05,
    "C\u03b1 RMSF, primer grip 225\u2013235 (\u00c5)": 0.05,
    "C\u03b1 RMSF, YMDD motif 183\u2013186 (\u00c5)": 0.05,
    "C\u03b1 RMSF, p66 fingers (\u00c5)": 0.05,
    "C\u03b1 RMSF, p66 thumb (\u00c5)": 0.05,
    "C\u03b1 RMSF, p66 connection (\u00c5)": 0.05,
}


def min_effect(column) -> float:
    """Minimum interpretable difference for one metric column."""
    if column in MIN_EFFECT:
        return MIN_EFFECT[column]
    unit = UNITS.get(column)
    if unit is None:
        raise KeyError(f"no unit known for {column!r}; add it to UNITS")
    return MIN_EFFECT_BY_UNIT[unit]


def resolved_mask(diff, diff_err, crit):
    """True where the shift clears BOTH the 95% criterion and its floor.

    `diff` and `diff_err` are frames whose columns are metric labels.
    """
    d = np.abs(diff.to_numpy(dtype=float))
    de = diff_err.reindex_like(diff).to_numpy(dtype=float)
    floor = np.array([min_effect(c) for c in diff.columns], dtype=float)
    return (d > crit * de) & (d > floor[None, :])

#: Category order of Table 1 in the manuscript.
CATEGORY_ORDER = ("Susceptible", "Resistant", "Uncertain")

#: Unit of each metric, used to group the columns. Putting an angle, a count
#: and a distance under one colour scale invites a comparison that cannot be
#: made: a +0.1 in Angstrom and a +0.1 in a contact count are different
#: quantities, and ordering them side by side implies otherwise. Columns are
#: therefore blocked by unit, with a rule and a header over each block.
UNITS = {
    "Lys103 C=O to DOR triazolinone N (\u00c5)": "\u00c5",
    "Lys103 side-chain polar atoms to DOR (\u00c5)": "\u00c5",
    "Val179 to DOR minimum distance (\u00c5)": "\u00c5",
    "Ser105 to DOR minimum distance (\u00c5)": "\u00c5",
    "Val106 to DOR minimum distance (\u00c5)": "\u00c5",
    "RT heavy-atom contacts < 4.0 \u00c5 of chlorocyanophenyl ring": "contacts",
    "Pyridinone moiety contacts < 4.0 \u00c5": "contacts",
    "Tyr188 / chlorocyanophenyl interplanar angle (\u00b0)": "degrees",
    "NNIBP volume (\u00c5\u00b3)": "\u00c5\u00b3",
    "Tyr181 side-chain \u03c71 (\u00b0)": "degrees",
    "Tyr188 side-chain \u03c71 (\u00b0)": "degrees",
    "C\u03b1 RMSF, NNIBP lining (\u00c5)": "\u00c5",
    "C\u03b1 RMSF, primer grip 225\u2013235 (\u00c5)": "\u00c5",
    "C\u03b1 RMSF, YMDD motif 183\u2013186 (\u00c5)": "\u00c5",
    "C\u03b1 RMSF, p66 fingers (\u00c5)": "\u00c5",
    "C\u03b1 RMSF, p66 thumb (\u00c5)": "\u00c5",
    "C\u03b1 RMSF, p66 connection (\u00c5)": "\u00c5",
    "\u03b29\u2013\u03b210 strand content (fraction)": "fraction",
}

#: Panel A columns, grouped by the DOR-RT interaction each one reports on,
#: in the order the first Results paragraph introduces them: the Lys103
#: hydrogen bond, the pyridinone packing on Val106, the chlorocyanophenyl
#: stacking on Tyr188, and then the three observables chosen from inspection
#: of the trajectories rather than from the crystal pose.
#:
#: Grouping by unit was the earlier scheme. It was defensible when colour was
#: a per-column ramp, but the fill is now categorical -- resolved above, below
#: or not at all -- so no cross-column magnitude comparison is on offer and the
#: unit blocks only separated quantities that belong together. A distance and
#: a contact count that report on the same interaction now sit side by side.
INTERACTION_GROUPS = {
    "Lys103 C=O to DOR triazolinone N (\u00c5)": "Lys103 H-bond",
    "Lys103 side-chain polar atoms to DOR (\u00c5)": "Lys103 H-bond",
    "Val106 to DOR minimum distance (\u00c5)": "Val106 packing",
    "Pyridinone moiety contacts < 4.0 \u00c5": "Val106 packing",
    "Ser105 to DOR minimum distance (\u00c5)": "Val106 packing",
    "Tyr188 / chlorocyanophenyl interplanar angle (\u00b0)": "Tyr188 stacking",
    "RT heavy-atom contacts < 4.0 \u00c5 of chlorocyanophenyl ring":
        "Tyr188 stacking",
    "Val179 to DOR minimum distance (\u00c5)": "Val179 contact",
}

GROUP_ORDER = ("Lys103 H-bond", "Val106 packing", "Tyr188 stacking",
               "Val179 contact")

#: Figure 5 columns, grouped by what each observable is meant to tell the
#: reader rather than by its unit. The chi1 dihedrals report on the rotamer
#: state of the two pocket-lining tyrosines; the volume and the lining RMSF
#: describe the pocket as a whole; the primer-grip RMSF is the one column
#: about RT away from the pocket.
RT_INTERPRETATION_GROUPS = {
    "Tyr181 side-chain \u03c71 (\u00b0)": "Local pocket remodeling",
    "Tyr188 side-chain \u03c71 (\u00b0)": "Local pocket remodeling",
    "NNIBP volume (\u00c5\u00b3)": "Global pocket properties",
    "C\u03b1 RMSF, NNIBP lining (\u00c5)": "Global pocket properties",
    "C\u03b1 RMSF, primer grip 225\u2013235 (\u00c5)":
        "RT conformational fluctuations",
}

RT_GROUP_ORDER = ("Local pocket remodeling", "Global pocket properties",
                  "RT conformational fluctuations")

#: Kept out of Figure 5. Neither the p66 fingers nor the p66 connection RMSF
#: is cited in the Results, and both are far enough from the pocket that a
#: resolved cell would not be interpretable on its own. They stay in the
#: panel C and D CSV outputs.
CD_EXCLUDE = ("C\u03b1 RMSF, p66 fingers (\u00c5)",
              "C\u03b1 RMSF, p66 connection (\u00c5)")


def order_by_rt_group(columns):
    """Figure 5 columns, blocked by intended interpretation."""
    cols, bounds, heads = [], [], []
    for g in RT_GROUP_ORDER:
        block = [c for c in RT_INTERPRETATION_GROUPS
                 if RT_INTERPRETATION_GROUPS[c] == g and c in columns]
        if not block:
            continue
        if cols:
            bounds.append(len(cols))
        heads.append((len(cols), len(cols) + len(block), g))
        cols += block
    cols += [c for c in columns if c not in cols]
    return cols, bounds, heads

#: Resolved in panel A, but kept out of the figure: the Val106-DOR minimum
#: distance moves over a 0.7 A range across the whole panel and clears its own
#: error in two genotypes by 0.27 and 0.28 A. The pyridinone contact count
#: reports the same interaction with more signal, so the column was costing a
#: reader's attention without adding a mechanism. It stays in the CSV outputs
#: and in Supplementary Table 4; only the drawn panel omits it.
#: The Val179-DOR minimum distance is likewise kept out of the drawn panel:
#: the same quantity is already shown in the Supplementary tables, so the
#: column is available to a reader who wants it without spending space in the
#: main figure. It too stays in the CSV outputs.
#: Val106 is no longer listed here: the minimum interpretable difference
#: removes it on its own, which is why that floor is not a post-hoc device.
PANEL_A_EXCLUDE = ("Val179 to DOR minimum distance (\u00c5)",)


def order_by_group(columns):
    """Columns grouped by the interaction they report on, with boundaries."""
    cols, bounds, heads = [], [], []
    for g in GROUP_ORDER:
        # Canonical order inside a block, so the columns appear in the order
        # the Results paragraph introduces them rather than in whatever order
        # the metric frame happens to carry.
        block = [c for c in INTERACTION_GROUPS
                 if INTERACTION_GROUPS[c] == g and c in columns]
        if not block:
            continue
        if cols:
            bounds.append(len(cols))
        heads.append((len(cols), len(cols) + len(block), g))
        cols += block
    cols += [c for c in columns if c not in cols]
    return cols, bounds, heads


#: Block order, widest-used unit first so the eye starts on the distances.
UNIT_ORDER = ("\u00c5", "\u00c5\u00b3", "degrees", "contacts", "fraction")

#: Resolved high, not resolved, resolved low. Three states, because three
#: replicates support three states: a continuous ramp invites the reader to
#: rank cells whose difference the data cannot resolve.
COLOUR_HIGH = "#C62828"
COLOUR_LOW = "#2166AC"
COLOUR_NONE = "#FFFFFF"


def order_by_unit(columns):
    """Columns grouped by unit, and the block boundaries between groups."""
    cols, bounds, heads = [], [], []
    for u in UNIT_ORDER:
        block = [c for c in columns if UNITS.get(c) == u]
        if not block:
            continue
        if cols:
            bounds.append(len(cols))
        heads.append((len(cols), len(cols) + len(block), u))
        cols += block
    cols += [c for c in columns if c not in cols]
    return cols, bounds, heads


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


def absolute_vs_wt(per_rep, state, metrics):
    """Each genotype's own value in one state, with WT's, and the shift.

    The cells print the measurement itself rather than its distance from wild
    type, which spares the reader doing the addition to recover a pocket
    volume or a bond length. The colour still answers the comparative
    question, so both are available at once: what the quantity is, and whether
    it is resolved against wild type.
    """
    d = per_rep[per_rep.structure_type == state]
    val, err, dif, dif_err = {}, {}, {}, {}
    for label, col in metrics.items():
        if col not in d or d[col].isna().all():
            continue
        st = _stats(d, col, label in CIRCULAR)
        if "WT" not in st.index or not np.isfinite(st.loc["WT", "mean"]):
            continue
        val[label] = st["mean"]
        err[label] = st["sem"]
        dif[label] = _wrap(st["mean"] - st.loc["WT", "mean"], label in CIRCULAR)
        dif_err[label] = np.sqrt(st["sem"] ** 2 + float(st.loc["WT", "sem"]) ** 2)
    return (pd.DataFrame(val), pd.DataFrame(err),
            pd.DataFrame(dif), pd.DataFrame(dif_err))


def induced_fit_absolute(per_rep, metrics):
    """Each genotype's own induced fit (holo - apo), and its shift from WT's."""
    val, err, dif, dif_err = {}, {}, {}, {}
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
        val[label] = fit
        err[label] = fit_sem
        dif[label] = _wrap(fit - fit["WT"], circ)
        dif_err[label] = np.sqrt(fit_sem ** 2 + float(fit_sem["WT"]) ** 2)
    return (pd.DataFrame(val), pd.DataFrame(err),
            pd.DataFrame(dif), pd.DataFrame(dif_err))


def induced_fit_vs_wt(per_rep, metrics):
    """How a genotype's induced fit differs from wild type's.

    The induced fit is delta = holo - apo: how far binding DOR moves a metric.
    Referencing that to wild type gives (holo - apo)_mut - (holo - apo)_WT, a
    difference of differences. Zero means the mutation shifts the bound and
    unliganded states by the same geometric amount. This does not establish
    cancellation of binding free energies. Four means are subtracted, so four
    SEMs add in quadrature.
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


def draw(ax, shown, shown_err, diff, diff_err, title, ylabels, show_y, crit,
         title_fs=44, title_pad=78):
    """Absolute values in the cells; colour says whether they beat wild type.

    `shown` carries each genotype's own measurement and `diff` its distance
    from wild type. The first row is wild type itself, which has no shift to
    test and is therefore never filled.

    `crit` is the number of standard errors the SHIFT must exceed for its 95%
    interval to exclude zero at this panel's degrees of freedom.
    """
    from matplotlib.colors import ListedColormap, BoundaryNorm

    v = shown.astype(float)
    ve = shown_err.reindex_like(v).astype(float)
    d = diff.reindex(v.index).reindex(columns=v.columns).astype(float)
    de = diff_err.reindex(v.index).reindex(columns=v.columns).astype(float)

    resolved = resolved_mask(d, de, crit)
    state = np.where(np.isfinite(d.to_numpy()) & resolved,
                     np.sign(d.to_numpy()), 0.0)

    cmap = ListedColormap([COLOUR_LOW, COLOUR_NONE, COLOUR_HIGH])
    norm = BoundaryNorm([-1.5, -0.5, 0.5, 1.5], cmap.N)
    im = ax.imshow(state, cmap=cmap, norm=norm, aspect="auto")

    # Decimals from the absolute values, which are larger than the shifts
    # were: a pocket volume of 228 does not want two decimal places.
    dec = {c: _decimals(float(np.nanmax(np.abs(v[c].to_numpy())))) for c in v.columns}
    for i in range(v.shape[0]):
        for j, c in enumerate(v.columns):
            x = v.iat[i, j]
            if not np.isfinite(x):
                ax.text(j, i, "\u00b7", ha="center", va="center", color="0.6",
                        fontsize=32)
                continue
            colour = "white" if state[i, j] != 0 else "0.12"
            ax.text(j, i - 0.15, f"{x:.{dec[c]}f}", ha="center", va="center",
                    fontsize=30, color=colour)
            y = ve.iat[i, j]
            if np.isfinite(y):
                ax.text(j, i + 0.24, f"\u00b1{y:.{dec[c]}f}", ha="center",
                        va="center", fontsize=23, color=colour, alpha=0.9)

    ax.set_xticks(range(v.shape[1]))
    ax.set_xticklabels([textwrap.fill(str(c), 26) for c in v.columns],
                       rotation=45, ha="right", rotation_mode="anchor",
                       fontsize=26)
    ax.set_yticks(range(v.shape[0]))
    ax.set_yticklabels(ylabels if show_y else [], fontsize=32)
    ax.set_title(title, fontsize=title_fs, pad=title_pad)
    ax.set_xticks(np.arange(-0.5, v.shape[1], 1), minor=True)
    ax.set_yticks(np.arange(-0.5, v.shape[0], 1), minor=True)
    ax.grid(which="minor", color="0.85", lw=1.0)
    ax.tick_params(which="minor", length=0)
    # Wild type is the reference, not a comparison: rule it off from the panel.
    ax.axhline(0.5, color="black", lw=3.0)
    return im


def draw_unit_blocks(ax, columns, n_rows, order=None, header_fs=30,
                     header_pad=12):
    """Rule between column blocks, and a header naming each block.

    ``order`` picks the blocking scheme: :func:`order_by_unit` by default, or
    :func:`order_by_group` for the panel that blocks by interaction instead.
    """
    _cols, bounds, heads = (order or order_by_unit)(list(columns))
    for b in bounds:
        ax.axvline(b - 0.5, color="black", lw=3.0)
    lines = 1
    for lo, hi, unit in heads:
        # A header has to fit over its own block, and a block can be one
        # column wide. Roughly 11 characters fit across a cell at this font,
        # so the text wraps to the block's width instead of running into its
        # neighbour -- which is exactly what the group names did at first.
        text = textwrap.fill(str(unit), max(9, int(11 * (hi - lo))),
                             break_long_words=False,
                             break_on_hyphens=False)
        lines = max(lines, text.count("\n") + 1)
        # Axes fraction 1.0 is the TOP of an imshow panel; 0 put these headers
        # underneath the last row of cells, where they overlapped the data.
        ax.annotate(text, xy=((lo + hi - 1) / 2.0, 1.0),
                    xycoords=("data", "axes fraction"),
                    xytext=(0, header_pad), textcoords="offset points",
                    ha="center", va="bottom", fontsize=header_fs,
                    fontweight="bold", color="0.25", clip_on=False,
                    linespacing=1.15, multialignment="center")
    # Height the headers occupy, in points, so the caller can stack the key
    # and the title above them instead of guessing a clearance.
    return header_pad + lines * header_fs * 1.2


def draw_key(ax, crit, pts_above=58, fontsize=25, inline=True):
    """The three-state colour key as a single row above the panel.

    It used to sit outside the right edge, which bought a column of white
    space the height of the whole panel for three lines of text. Placed over
    the panel instead, under the title, it costs one row and the figure loses
    several inches of width.

    The offset is in POINTS from the top of the axes, not a fraction of it:
    these panels are over twenty inches tall, so an axes fraction that looks
    right on one panel throws the key into the title on another.
    """
    from matplotlib.patches import Patch
    from matplotlib.transforms import ScaledTranslation

    handles = []
    if inline:
        # An invisible patch renders its label as leading text, which keeps
        # the criterion on the same line instead of stacking a legend title
        # above the row.
        handles.append(Patch(facecolor="none", edgecolor="none",
                             label=f"Resolved vs WT ({crit:.3g}\u00d7 SEM and "
                                   f"minimum effect):"))
    handles += [
        Patch(facecolor=COLOUR_HIGH, edgecolor="0.4", label="above WT"),
        Patch(facecolor=COLOUR_LOW, edgecolor="0.4", label="below WT"),
        Patch(facecolor=COLOUR_NONE, edgecolor="0.4",
              label="not resolved")]
    trans = ax.transAxes + ScaledTranslation(
        0.0, pts_above / 72.0, ax.figure.dpi_scale_trans)
    return ax.legend(handles=handles, loc="lower center",
                     bbox_to_anchor=(0.5, 1.0), bbox_transform=trans,
                     ncol=len(handles), frameon=False, fontsize=fontsize,
                     handlelength=1.5, handletextpad=0.5,
                     columnspacing=1.4, borderpad=0.0, borderaxespad=0.0)


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--per-replicate", type=Path,
                   default=Path("results/analysis/structural_metrics/tables/"
                                "structural_metrics_per_replicate.csv"))
    p.add_argument("--table3-csv", type=Path,
                   default=None,
                   help="Optional legacy category CSV; defaults to the shared panel definition.")
    p.add_argument("--manifest", type=Path, default=Path("manifests/md_manifest.csv"))
    p.add_argument("--output-dir", type=Path,
                   default=Path("results/analysis/structural_metrics"))
    args = p.parse_args(argv)

    plt = _style()
    from matplotlib.transforms import blended_transform_factory
    per_rep = pd.read_csv(args.per_replicate)
    man = pd.read_csv(args.manifest)[["mutation", "fold_reduction"]].drop_duplicates()
    fold = dict(zip(man.mutation, man.fold_reduction))
    from nnrti.cli.build_supplementary_table_4 import CATEGORY
    cat = dict(CATEGORY) if args.table3_csv is None else dict(
        pd.read_csv(args.table3_csv)[["Genotype", "Category"]].values)


    cat = {g: c for g, c in cat.items() if g == "WT" or c in CATEGORY_ORDER}

    dor_a, dor_ae, dor_v, dor_e = absolute_vs_wt(per_rep, "holo", DOR_METRICS)
    rt_a, rt_ae, rt_v, rt_e = induced_fit_absolute(per_rep, RT_METRICS)
    # RT conformation referenced to wild type WITHIN each state. This answers
    # a different question from the induced-fit panel: that one asks how
    # differently the mutation lets doravirine reshape the enzyme, these ask
    # how the mutant enzyme differs from wild type in a single state, with no
    # cancellation between the two. A genotype can be flat in the induced-fit
    # panel and far from wild type in both of these.
    rtapo_a, rtapo_ae, rtapo_v, rtapo_e = absolute_vs_wt(per_rep, "apo", RT_METRICS)
    rtholo_a, rtholo_ae, rtholo_v, rtholo_e = absolute_vs_wt(per_rep, "holo",
                                                             RT_METRICS)

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
    # The fold reduction still ORDERS the rows within a resistance category,
    # but is no longer printed beside each genotype: it is in Table 1 and in
    # Figure 2, and repeating it here only crowded the axis.
    labels = list(mutants)

    def rule_and_label(ax, rows, offset=0):
        """Category rules and side labels. `offset` skips the wild-type row."""
        cats = [cat[g] for g in rows]
        for i in range(1, len(cats)):
            if cats[i] != cats[i - 1]:
                ax.axhline(i + offset - 0.5, color="black", lw=2.4)
        seen, start_i = cats[0], 0
        for i in range(1, len(cats) + 1):
            if i == len(cats) or cats[i] != seen:
                # Offset in POINTS, not axes fraction: a fraction of the
                # axes width is ~5 inches on a 20-inch-wide panel, which
                # stranded these labels far out in the margin.
                ax.annotate(seen, xy=(0, (start_i + i - 1) / 2.0 + offset),
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
        sig = pd.DataFrame(resolved_mask(v, e, crit),
                           index=v.index, columns=v.columns).sum()
        drop = [c for c in v.columns if sig[c] == 0]
        if drop:
            print(f"{tag}: dropped, no genotype clears both the 95% "
                  f"criterion ({crit:.3g}x SEM) and the minimum "
                  f"interpretable difference -- {', '.join(drop)}")
        keep = [c for c in v.columns if sig[c] > 0]
        return val[keep], err[keep]

    # The CSV outputs are the data record and must not lose a metric just
    # because the drawn panel did: a column dropped for showing nothing is
    # itself a result, and a reader checking it needs the numbers.
    unfiltered = {"A": (dor_v, dor_e, dor_a, dor_ae),
                  "B": (rt_v, rt_e, rt_a, rt_ae),
                  "C": (rtapo_v, rtapo_e, rtapo_a, rtapo_ae),
                  "D": (rtholo_v, rtholo_e, rtholo_a, rtholo_ae)}

    dor_v, dor_e = keep_resolved(dor_v, dor_e, mutants, "A")
    rt_v, rt_e = keep_resolved(rt_v, rt_e, mutants, "B")
    rtapo_v, rtapo_e = keep_resolved(rtapo_v, rtapo_e, mutants, "C")
    rtholo_v, rtholo_e = keep_resolved(rtholo_v, rtholo_e, mutants, "D")

    # Titles no longer say "shift from wild type": the cells carry the
    # measurement itself, and only the fill is comparative.
    panels = [
        ("A", dor_a, dor_ae, dor_v, dor_e,
         "RT\u2013DOR specific metrics (bound)",
         "figure3A_rt_dor_contacts", "figure3A_dor"),
        ("B", rt_a, rt_ae, rt_v, rt_e,
         "RT-specific metrics: induced fit (bound \u2212 unliganded)",
         "figure3B_rt_induced_fit", "figure3B_rt"),
        ("C", rtapo_a, rtapo_ae, rtapo_v, rtapo_e,
         "RT-specific metrics (unliganded)",
         "figure3C_rt_apo_vs_wt", "figure3C_rt_apo"),
        ("D", rtholo_a, rtholo_ae, rtholo_v, rtholo_e,
         "RT-specific metrics (bound)",
         "figure3D_rt_holo_vs_wt", "figure3D_rt_holo"),
    ]
    written = []
    for tag, absv, abse, val, err, title, stem, tstem in panels:
        if val.shape[1] == 0:
            print(f"{tag}: no metric has a resolved cell; panel not drawn")
            continue
        # Panel B is Supplementary Figure 7 and is read against Figure 5, so
        # it drops the same two RMSF columns and blocks by the same
        # interpretation groups.
        drop = {"A": PANEL_A_EXCLUDE, "B": CD_EXCLUDE}.get(tag, ())
        panel_cols = [c for c in val.columns if c not in drop]
        blocking = {"A": order_by_group,
                    "B": order_by_rt_group}.get(tag, order_by_unit)
        cols, _b, _h = blocking(panel_cols)
        # Wild type heads the panel as the reference row.
        rows = ["WT"] + mutants
        v = val.reindex(rows)[cols]
        e = err.reindex(rows)[cols]
        av = absv.reindex(rows)[cols]
        ae = abse.reindex(rows)[cols]
        # The CSVs keep every resolved metric, including one the drawn panel
        # leaves out, so the figure's omission never becomes a hole in the
        # data record.
        u_val, u_err, u_absv, u_abse = unfiltered[tag]
        all_cols, _b2, _h2 = blocking(list(u_val.columns))
        v_all = u_val.reindex(rows)[all_cols]
        e_all = u_err.reindex(rows)[all_cols]
        av_all = u_absv.reindex(rows)[all_cols]
        ae_all = u_abse.reindex(rows)[all_cols]
        # Canvas deliberately NOT scaled with the font. Growing both together
        # leaves the type the same size relative to the cells and only makes
        # the file enormous -- doubling the fonts and the figure gave a 46 inch
        # wide image that read exactly as before.
        # Shrinking the canvas while the fonts stay fixed in points makes the
        # type larger relative to the cells, which is the goal, and keeps the
        # render under 8000 px on its long side so it can be viewed anywhere.
        fig, ax = plt.subplots(figsize=(2.05 * v.shape[1] + 5.0,
                                        0.92 * len(rows)
                                        + (6.2 if tag == "A" else 5.0)))
        im = draw(ax, av, ae, v, e, title, ["WT"] + labels, True,
                  SIGNIFICANCE_BY_PANEL[tag])
        rule_and_label(ax, mutants, offset=1)
        header_fs = 28 if tag == "A" else 30
        hdr_pts = draw_unit_blocks(ax, v.columns, len(mutants),
                                   order=blocking, header_fs=header_fs)
        # A three-state key, not a colour bar: the fill says only whether the
        # 95% interval excludes zero and in which direction. Magnitudes are in
        # the cells, in each metric's own units, which is the only place they
        # can be compared honestly.
        crit = SIGNIFICANCE_BY_PANEL[tag]
        if tag in ("A", "B"):
            # Title at the top, then the key on one line, then the group
            # headers sitting directly on the cells they label.
            key_pts = hdr_pts + 26
            title_fs = 44 if tag == "A" else 38
            draw_key(ax, crit, pts_above=key_pts)
            ax.set_title(title, fontsize=title_fs,
                         pad=key_pts + 28 * 1.9 + 26)
        else:
            from matplotlib.patches import Patch
            ax.legend(handles=[
                Patch(facecolor=COLOUR_HIGH, edgecolor="0.4", label="above WT"),
                Patch(facecolor=COLOUR_LOW, edgecolor="0.4", label="below WT"),
                Patch(facecolor=COLOUR_NONE, edgecolor="0.4",
                      label="not resolved")],
                loc="upper left", bbox_to_anchor=(1.01, 1.0), frameon=False,
                fontsize=28, title=f"Resolved vs WT\n({crit:.3g}\u00d7 SEM and minimum effect)",
                title_fontsize=28)
        f = out / "plots" / f"{stem}.png"
        # 260, not 300: panel A is 27 inches wide and renders past the
        # 8000 px that browsers and phone clients will display.
        fig.savefig(f, dpi=260, bbox_inches="tight", pad_inches=0.25)
        plt.close(fig)
        written.append(f)
        v_all.loc[mutants].round(3).to_csv(
            out / "tables" / f"{tstem}_shift_vs_wt.csv")
        e_all.loc[mutants].round(3).to_csv(
            out / "tables" / f"{tstem}_shift_sem.csv")
        av_all.round(3).to_csv(out / "tables" / f"{tstem}_absolute.csv")
        ae_all.round(3).to_csv(out / "tables" / f"{tstem}_absolute_sem.csv")
        vm, em = v.loc[mutants], e.loc[mutants]
        n = int(np.isfinite(vm.to_numpy()).sum())
        sig = int(resolved_mask(vm, em, SIGNIFICANCE_BY_PANEL[tag]).sum())
        print(f"{tag}: {vm.shape[0]} genotypes x {vm.shape[1]} metrics, "
              f"{sig}/{n} cells resolved against WT")


    # ---- the two single-state panels, side by side --------------------------
    # One figure rather than two: the rows are the same genotypes in the same
    # order, so the y axis, the category rules and the colour key are each
    # drawn once, and the reader compares a metric between the two states by
    # moving sideways instead of between pages. Columns are the UNION of what
    # the two panels resolve, so a metric kept in one state still appears in
    # the other -- all white there, which is itself the answer.
    cd_cols = [c for c in RT_METRICS
               if c in set(rtholo_v.columns) | set(rtapo_v.columns)
               and c not in CD_EXCLUDE]
    cd_cols, _b, _h = order_by_rt_group(cd_cols)
    rows = ["WT"] + mutants
    # Narrower cells than the single panels use, and a lower dpi, because
    # twice the columns at the single-panel scale renders past 10,000 px and
    # will not display on a phone or in a browser.
    # NOT sharey: shared axes share one set of tick labels, so blanking the
    # right panel's y axis blanked the left panel's genotypes too. The two
    # panels have the same number of rows and therefore the same limits
    # anyway.
    fig, axes = plt.subplots(1, 2,
                             figsize=(1.95 * 2 * len(cd_cols) + 9.0,
                                      0.92 * len(rows) + 6.0),
                             gridspec_kw={"wspace": 0.05})
    for k, (absv, abse, val, err, sub_title, tag) in enumerate((
            (rtholo_a, rtholo_ae, rtholo_v, rtholo_e, "Bound (RT+DOR)", "D"),
            (rtapo_a, rtapo_ae, rtapo_v, rtapo_e, "Unliganded (RT only)", "C"))):
        ax = axes[k]
        av = absv.reindex(rows).reindex(columns=cd_cols)
        ae = abse.reindex(rows).reindex(columns=cd_cols)
        v = val.reindex(rows).reindex(columns=cd_cols)
        e = err.reindex(rows).reindex(columns=cd_cols)
        draw(ax, av, ae, v, e, sub_title, ["WT"] + labels, k == 0,
             SIGNIFICANCE_BY_PANEL[tag], title_fs=38, title_pad=64)
        hdr_pts = draw_unit_blocks(ax, cd_cols, len(mutants),
                                   order=order_by_rt_group,
                                   header_fs=26)
        # Re-set the subtitle now that the headers have been measured: at
        # draw()'s fixed pad it landed on top of them.
        ax.set_title(sub_title, fontsize=38, pad=hdr_pts + 22)
        if k == 0:
            rule_and_label(ax, mutants, offset=1)
        else:
            cats = [cat[g] for g in mutants]
            for i in range(1, len(cats)):
                if cats[i] != cats[i - 1]:
                    ax.axhline(i + 0.5, color="black", lw=2.4)
    # The key used to sit outside the right panel, which widened the figure
    # by its own width for three swatches. On one centred line under the
    # title it costs a row of text instead of a column.
    #
    # Everything above the cells is stacked from the measured header height,
    # in points, and centred on the figure in x: placing any of it at a fixed
    # figure fraction puts it on top of the panel subtitles, whose own height
    # depends on how many lines the group headers wrapped to.
    from matplotlib.patches import Patch
    from matplotlib.transforms import ScaledTranslation, blended_transform_factory
    crit = SIGNIFICANCE_BY_PANEL["D"]
    handles = [
        Patch(facecolor="none", edgecolor="none",
              label=f"Resolved vs WT ({crit:.3g}\u00d7 SEM and "
                    f"minimum effect):"),
        Patch(facecolor=COLOUR_HIGH, edgecolor="0.4", label="above WT"),
        Patch(facecolor=COLOUR_LOW, edgecolor="0.4", label="below WT"),
        Patch(facecolor=COLOUR_NONE, edgecolor="0.4", label="not resolved")]

    def _above(pts):
        """Figure-centred x, `pts` above the top of the panels in y."""
        return blended_transform_factory(
            fig.transFigure,
            axes[0].transAxes + ScaledTranslation(
                0.0, pts / 72.0, fig.dpi_scale_trans))

    sub_pts = hdr_pts + 22 + 38 * 1.25          # top of the panel subtitles
    key_pts = sub_pts + 22
    fig.legend(handles=handles, loc="lower center",
               bbox_to_anchor=(0.5, 1.0), bbox_transform=_above(key_pts),
               frameon=False, fontsize=26, ncol=len(handles),
               handlelength=1.5, handletextpad=0.5, columnspacing=1.4,
               borderpad=0.0, borderaxespad=0.0)
    fig.text(0.5, 1.0, "RT-specific metrics", fontsize=48, ha="center",
             va="bottom", transform=_above(key_pts + 26 * 1.3 + 24))
    f = out / "plots" / "figure3CD_rt_by_state.png"
    fig.savefig(f, dpi=185, bbox_inches="tight", pad_inches=0.3)
    plt.close(fig)
    written.append(f)
    print(f"CD: combined sheet, {len(cd_cols)} metrics x {len(rows)} rows")

    for f in written:
        print(f"wrote {f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
