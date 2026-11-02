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
from pathlib import Path

import numpy as np
import pandas as pd

#: Panel A: contacts between RT and DOR. Undefined without a ligand.
DOR_METRICS = {
    "d(103 C=O···N4x)": "res103_CO_to_triazolinone_N_A",
    "d(103 sc···DOR)": "res103_sidechain_to_DOR_A",
    "N(chlorophenyl)": "chlorocyanophenyl_contacts_n",
    "θ(Y188/chl)": "Tyr188_chlorophenyl_angle_deg",
    "d(179···DOR)": "res179_to_DOR_A",
    "d(S105···DOR)": "Ser105_to_DOR_A",
    "d(106···DOR)": "res106_to_DOR_A",
    "N(pyridinone)": "pyridinone_contacts_n",
}

#: Panel B: RT's own conformation, measurable with or without the ligand.
RT_METRICS = {
    "V(NNIBP)": "nnibp_volume_A3",
    "χ1 Tyr181": "chi1_Tyr181_deg",
    "χ1 Tyr188": "chi1_Tyr188_deg",
    "RMSF pocket": "rmsf_A_nnibp_lining",
    "RMSF primer grip": "rmsf_A_primer_grip",
    "RMSF YMDD": "rmsf_A_ymdd",
    "RMSF fingers": "rmsf_A_fingers",
    "RMSF thumb": "rmsf_A_thumb",
    "RMSF connection": "rmsf_A_connection",
    "strand β9–β10": "strand_frac_beta9_beta10",
}

#: Dihedrals: mean and difference must both be taken on the circle.
CIRCULAR = {"χ1 Tyr181", "χ1 Tyr188"}

UNITS = {
    "d(103 C=O···N4x)": "Å", "d(103 sc···DOR)": "Å", "N(chlorophenyl)": "count",
    "θ(Y188/chl)": "°", "d(179···DOR)": "Å", "d(S105···DOR)": "Å",
    "d(106···DOR)": "Å", "N(pyridinone)": "count", "V(NNIBP)": "Å³",
    "χ1 Tyr181": "°", "χ1 Tyr188": "°", "RMSF pocket": "Å",
    "RMSF primer grip": "Å", "RMSF YMDD": "Å", "RMSF fingers": "Å",
    "RMSF thumb": "Å", "RMSF connection": "Å", "strand β9–β10": "frac",
}

#: Category order of Table 1 in the manuscript.
CATEGORY_ORDER = ("Susceptible", "Resistant", "Uncertain")


def _style():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({
        "font.size": 11, "axes.labelsize": 11, "axes.titlesize": 13,
        "xtick.labelsize": 10, "ytick.labelsize": 10,
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


def holo_minus_apo(per_rep, metrics):
    """Induced fit: how far binding DOR moves each metric, per genotype."""
    val, err = {}, {}
    for label, col in metrics.items():
        if col not in per_rep or per_rep[col].isna().all():
            continue
        h = _stats(per_rep[per_rep.structure_type == "holo"], col, label in CIRCULAR)
        a = _stats(per_rep[per_rep.structure_type == "apo"], col, label in CIRCULAR)
        common = h.index.intersection(a.index)
        val[label] = _wrap(h.loc[common, "mean"] - a.loc[common, "mean"],
                           label in CIRCULAR)
        err[label] = np.sqrt(h.loc[common, "sem"] ** 2 + a.loc[common, "sem"] ** 2)
    return pd.DataFrame(val), pd.DataFrame(err)


def _decimals(largest: float) -> int:
    return 0 if largest >= 10 else (1 if largest >= 1 else 2)


def draw(ax, val, err, title, ylabels, show_y):
    scaled = val.astype(float).copy()
    for c in scaled.columns:
        m = float(np.nanmax(np.abs(scaled[c].to_numpy())))
        scaled[c] = scaled[c] / m if m else 0.0
    ax.imshow(scaled.to_numpy(dtype=float), cmap="RdBu_r", vmin=-1, vmax=1,
              aspect="auto")

    dec = {c: _decimals(float(np.nanmax(np.abs(val[c].to_numpy()))))
           for c in val.columns}
    for i in range(val.shape[0]):
        for j, c in enumerate(val.columns):
            v = val.iat[i, j]
            if not np.isfinite(v):
                ax.text(j, i, "·", ha="center", va="center", color="0.6", fontsize=10)
                continue
            colour = "white" if abs(float(scaled.iat[i, j])) > 0.65 else "0.12"
            ax.text(j, i - 0.16, f"{v:+.{dec[c]}f}", ha="center", va="center",
                    fontsize=8.5, color=colour)
            e = err.iat[i, j]
            if np.isfinite(e):
                ax.text(j, i + 0.22, f"±{e:.{dec[c]}f}", ha="center", va="center",
                        fontsize=6.8, color=colour, alpha=0.9)

    ax.set_xticks(range(val.shape[1]))
    ax.set_xticklabels([f"{c}  ({UNITS.get(c, '')})" for c in val.columns],
                       rotation=45, ha="right", rotation_mode="anchor", fontsize=10)
    ax.set_yticks(range(val.shape[0]))
    ax.set_yticklabels(ylabels if show_y else [], fontsize=10)
    ax.set_title(title, fontsize=13, pad=12)
    ax.set_xticks(np.arange(-0.5, val.shape[1], 1), minor=True)
    ax.set_yticks(np.arange(-0.5, val.shape[0], 1), minor=True)
    ax.grid(which="minor", color="white", lw=1.0)
    ax.tick_params(which="minor", length=0)


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

    # A genotype that is simulated but absent from Table 1 has no category, and
    # silently sorting it into one would misrepresent the manuscript. Name it.
    missing = sorted({g for g in per_rep.genotype.unique()
                      if g not in cat and g != "WT"})
    if missing:
        print(f"NOT in Table 1, so shown in a separate block: {', '.join(missing)}")

    dor_v, dor_e = shift_vs_wt(per_rep, "holo", DOR_METRICS)
    rt_v, rt_e = holo_minus_apo(per_rep, RT_METRICS)

    # Category first, then fold-change within category, as Table 1 groups them.
    def category(g: str) -> str:
        return cat.get(g, "Not in Table 1")

    def key(g):
        c = category(g)
        rank = CATEGORY_ORDER.index(c) if c in CATEGORY_ORDER else len(CATEGORY_ORDER)
        return (rank, np.nan_to_num(fold.get(g, np.nan), nan=np.inf))

    mutants = sorted((g for g in dor_v.index if g != "WT"), key=key)
    labels_a = [f"{g}   {fold.get(g):g}×" if np.isfinite(fold.get(g, np.nan)) else g
                for g in mutants]
    # Panel B keeps WT: the question is whether a genotype's induced fit
    # DIFFERS from wild type's, which needs wild type on the page.
    rt_order = ["WT"] + [g for g in mutants if g in rt_v.index]
    labels_b = ["WT (reference)"] + [l for g, l in zip(mutants, labels_a)
                                     if g in rt_v.index]

    fig, axes = plt.subplots(
        1, 2, figsize=(26, 0.62 * len(mutants) + 6.5),
        gridspec_kw={"width_ratios": [dor_v.shape[1], rt_v.shape[1]], "wspace": 0.30})
    draw(axes[0], dor_v.loc[mutants], dor_e.reindex(mutants),
         "A   RT–DOR contacts (holo), shift from wild type", labels_a, True)
    draw(axes[1], rt_v.loc[rt_order], rt_e.reindex(rt_order),
         "B   RT conformation, Δ = holo − apo (induced fit)", labels_b, True)

    # Rule between resistance categories, and name each block in the margin.
    def rule_and_label(ax, rows, first_is_wt, label=True):
        from matplotlib.transforms import blended_transform_factory
        cats = ["WT"] + [category(g) for g in rows[1:]] if first_is_wt \
            else [category(g) for g in rows]
        for i in range(1, len(cats)):
            if cats[i] != cats[i - 1]:
                ax.axhline(i - 0.5, color="black", lw=2.4)
        seen, start = cats[0], 0
        for i in range(1, len(cats) + 1):
            if i == len(cats) or cats[i] != seen:
                if seen != "WT" and label:
                    # x in axes fraction so the label clears the genotype tick
                    # labels whatever their length; y in data so it centres on
                    # the block. A single transform cannot do both.
                    tr = blended_transform_factory(ax.transAxes, ax.transData)
                    ax.text(-0.30, (start + i - 1) / 2.0, seen, rotation=90,
                            ha="center", va="center", fontsize=12,
                            fontweight="bold", color="0.25",
                            transform=tr, clip_on=False)
                if i < len(cats):
                    seen, start = cats[i], i
    # Only the left panel is labelled: both share the row order, and panel B's
    # labels would sit in the gap and read as part of panel A.
    rule_and_label(axes[0], mutants, False, label=True)
    rule_and_label(axes[1], rt_order, True, label=False)

    out = args.output_dir
    (out / "plots").mkdir(parents=True, exist_ok=True)
    (out / "tables").mkdir(parents=True, exist_ok=True)
    f = out / "plots" / "figure3_structural_heatmap.png"
    fig.savefig(f, dpi=300)
    plt.close(fig)

    dor_v.loc[mutants].round(3).to_csv(out / "tables" / "figure3A_dor_shift_vs_wt.csv")
    dor_e.reindex(mutants).round(3).to_csv(out / "tables" / "figure3A_dor_shift_sem.csv")
    rt_v.loc[rt_order].round(3).to_csv(out / "tables" / "figure3B_rt_holo_minus_apo.csv")
    rt_e.reindex(rt_order).round(3).to_csv(out / "tables" / "figure3B_rt_holo_minus_apo_sem.csv")

    sig_a = int((dor_v.loc[mutants].abs() > 2 * dor_e.reindex(mutants)).to_numpy().sum())
    n_a = int(np.isfinite(dor_v.loc[mutants].to_numpy()).sum())
    body = rt_v.loc[mutants]
    sig_b = int((body.abs() > 2 * rt_e.reindex(mutants)).to_numpy().sum())
    n_b = int(np.isfinite(body.to_numpy()).sum())
    print(f"A  RT–DOR contacts : {sig_a}/{n_a} cells exceed 2x SEM")
    print(f"B  RT induced fit  : {sig_b}/{n_b} mutant cells exceed 2x SEM")
    print(f"\nwrote {f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
