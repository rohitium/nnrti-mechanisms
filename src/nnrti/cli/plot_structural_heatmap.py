#!/usr/bin/env python3
"""Figure 3: WT-referenced structural metrics across the panel.

Reviewer 1 (comment 4) asked for a WT-referenced heatmap summarising the
structural metrics, because they are currently spread across dense tables and
patterns in DOR contacts, pocket size and remodelling are hard to see.

Each cell is a genotype's shift from wild type IN ITS OWN UNITS, printed in
the cell. Colour is scaled within each column, to that metric's largest shift
across the panel, so a column can be read down at a glance; it carries no
meaning across columns, and the number is the quantity.

An earlier version divided each shift by the pooled replicate standard
deviation, to put every metric on one scale. That inverted the physics. G190E's
pocket expansion -- the largest structural change in the panel at +55.8 A^3,
and the one the reviewer singles out -- came out at +3.4 SD and rendered pale,
while a 2.84 A shift in d(179...DOR) came out at +10.6 SD and rendered dark,
purely because contact distances have tighter replicate scatter than pocket
volumes. An effect size answers "is this bigger than the noise"; this figure
has to answer "how much did the structure move", and those are different
questions.

Rows are ordered by DOR fold resistance, so a reader can see whether
structural change tracks phenotype at all.

Usage
-----
    PYTHONPATH=src python -m nnrti.cli.plot_structural_heatmap
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

#: Display name -> column, in the order they appear across the heatmap.
HOLO_METRICS = {
    "d(103 C=O···N4x)": "res103_CO_to_triazolinone_N_A",
    "d(103 sc···DOR)": "res103_sidechain_to_DOR_A",
    "N(chlorophenyl)": "chlorocyanophenyl_contacts_n",
    "θ(Y188/chl)": "Tyr188_chlorophenyl_angle_deg",
    "d(179···DOR)": "res179_to_DOR_A",
    "d(S105···DOR)": "Ser105_to_DOR_A",
    "d(106···DOR)": "res106_to_DOR_A",
    "N(pyridinone)": "pyridinone_contacts_n",
}
POCKET_METRICS = {
    "V(NNIBP)": "nnibp_volume_A3",
    "χ1 Tyr181": "chi1_Tyr181_deg",
    "χ1 Tyr188": "chi1_Tyr188_deg",
}
PROTEIN_METRICS = {
    "RMSF pocket": "rmsf_A_nnibp_lining",
    "RMSF primer grip": "rmsf_A_primer_grip",
    "RMSF YMDD": "rmsf_A_ymdd",
    "RMSF fingers": "rmsf_A_fingers",
    "RMSF thumb": "rmsf_A_thumb",
    "RMSF connection": "rmsf_A_connection",
    "helix YMDD": "helix_frac_ymdd",
    "strand β9–β10": "strand_frac_beta9_beta10",
}


#: Metrics that are dihedral angles. Their panel mean and their difference
#: from WT must both be taken on the circle: a linear difference reported
#: Y318F's apo Tyr188 chi1 as +278 degrees, when the actual rotation is -82.
CIRCULAR = {"\u03c71 Tyr181", "\u03c71 Tyr188"}

#: Unit shown under each column header. The cell values are in these units.
UNITS = {
    "d(103 C=O\u00b7\u00b7\u00b7N4x)": "\u00c5", "d(103 sc\u00b7\u00b7\u00b7DOR)": "\u00c5",
    "N(chlorophenyl)": "count", "\u03b8(Y188/chl)": "\u00b0",
    "d(179\u00b7\u00b7\u00b7DOR)": "\u00c5", "d(S105\u00b7\u00b7\u00b7DOR)": "\u00c5",
    "d(106\u00b7\u00b7\u00b7DOR)": "\u00c5", "N(pyridinone)": "count",
    "V(NNIBP)": "\u00c5\u00b3", "\u03c71 Tyr181": "\u00b0", "\u03c71 Tyr188": "\u00b0",
    "RMSF pocket": "\u00c5", "RMSF primer grip": "\u00c5", "RMSF YMDD": "\u00c5",
    "RMSF fingers": "\u00c5", "RMSF thumb": "\u00c5", "RMSF connection": "\u00c5",
    "helix YMDD": "frac", "strand \u03b29\u2013\u03b210": "frac",
}



def _style():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({
        "font.size": 8, "axes.labelsize": 8, "axes.titlesize": 9,
        "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
        "figure.dpi": 300, "savefig.bbox": "tight",
    })
    return plt


def wt_referenced(per_rep: pd.DataFrame, state: str,
                  metrics: dict[str, str]) -> pd.DataFrame:
    """Shift from WT, in each metric's own units."""
    d = per_rep[per_rep.structure_type == state]
    raw = {}
    for label, col in metrics.items():
        if col not in d or d[col].isna().all():
            continue
        if label in CIRCULAR:
            rad = d.assign(_r=np.deg2rad(d[col]))
            means = rad.groupby("genotype").apply(
                lambda g: np.rad2deg(np.arctan2(np.sin(g._r).mean(),
                                                np.cos(g._r).mean())),
                include_groups=False)
        else:
            means = d.groupby("genotype")[col].mean()
        wt = means.get("WT")
        if wt is None or not np.isfinite(wt):
            continue
        if label in CIRCULAR:
            # Shortest rotation, so a dihedral that moved 82 degrees is never
            # reported as having moved 278 the other way round.
            delta = ((means - wt + 180.0) % 360.0) - 180.0
        else:
            delta = means - wt
        # A metric that never moves has nothing to show and no scale to show it
        # on; drop it rather than plot an all-zero column.
        if not np.isfinite(delta.drop(index="WT", errors="ignore")).any():
            continue
        if float(np.nanmax(np.abs(delta.drop(index="WT", errors="ignore")))) == 0.0:
            continue
        raw[label] = delta
    return pd.DataFrame(raw)


def _decimals(largest: float) -> int:
    """Enough precision to distinguish cells, no more."""
    if largest >= 10:
        return 0
    if largest >= 1:
        return 1
    return 2


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--per-replicate", type=Path,
                   default=Path("results/analysis/structural_metrics/tables/"
                                "structural_metrics_per_replicate.csv"))
    p.add_argument("--manifest", type=Path, default=Path("manifests/md_manifest.csv"))
    p.add_argument("--output-dir", type=Path,
                   default=Path("results/analysis/structural_metrics"))
    args = p.parse_args(argv)

    plt = _style()
    from matplotlib.colors import TwoSlopeNorm

    per_rep = pd.read_csv(args.per_replicate)
    man = pd.read_csv(args.manifest)[["mutation", "fold_reduction"]].drop_duplicates()
    fold = dict(zip(man.mutation, man.fold_reduction))

    all_metrics = {**HOLO_METRICS, **POCKET_METRICS, **PROTEIN_METRICS}
    raw_holo = wt_referenced(per_rep, "holo", all_metrics)
    apo_metrics = {**POCKET_METRICS, **PROTEIN_METRICS}
    raw_apo = wt_referenced(per_rep, "apo", apo_metrics)

    order = sorted((g for g in raw_holo.index if g != "WT"),
                   key=lambda g: (np.nan_to_num(fold.get(g, np.nan), nan=-1)))
    raw_holo = raw_holo.loc[order]
    raw_apo = raw_apo.reindex(order)

    panels = [("DOR-bound (holo)", raw_holo), ("Unliganded (apo)", raw_apo)]
    widths = [p[1].shape[1] for p in panels]
    fig, axes = plt.subplots(
        1, 2, figsize=(0.62 * sum(widths) + 3.4, 0.34 * len(order) + 2.4),
        gridspec_kw={"width_ratios": widths, "wspace": 0.05})

    for ax, (title, raw) in zip(axes, panels):
        # Colour is scaled per column: each metric spans its own largest shift,
        # so the gradient compares genotypes WITHIN a metric and never across.
        scaled = raw.copy().astype(float)
        for c in scaled.columns:
            m = float(np.nanmax(np.abs(scaled[c].to_numpy())))
            scaled[c] = scaled[c] / m if m else 0.0
        ax.imshow(scaled.to_numpy(dtype=float), cmap="RdBu_r",
                  vmin=-1, vmax=1, aspect="auto")

        dec = {c: _decimals(float(np.nanmax(np.abs(raw[c].to_numpy()))))
               for c in raw.columns}
        for i in range(raw.shape[0]):
            for j, c in enumerate(raw.columns):
                v = raw.iat[i, j]
                if not np.isfinite(v):
                    ax.text(j, i, "\u00b7", ha="center", va="center",
                            color="0.6", fontsize=7)
                    continue
                shade = abs(float(scaled.iat[i, j]))
                ax.text(j, i, f"{v:+.{dec[c]}f}", ha="center", va="center",
                        fontsize=5.5,
                        color="white" if shade > 0.65 else "0.15")

        ax.set_xticks(range(raw.shape[1]))
        ax.set_xticklabels([f"{c}\n{UNITS.get(c, '')}" for c in raw.columns],
                           rotation=90, ha="center", fontsize=6)
        ax.set_yticks(range(len(order)))
        ax.set_yticklabels([f"{g}  ({fold.get(g):g}\u00d7)" if np.isfinite(
            fold.get(g, np.nan)) else g for g in order] if ax is axes[0] else [])
        ax.set_title(title, fontsize=9)
        ax.set_xticks(np.arange(-0.5, raw.shape[1], 1), minor=True)
        ax.set_yticks(np.arange(-0.5, len(order), 1), minor=True)
        ax.grid(which="minor", color="white", lw=0.8)
        ax.tick_params(which="minor", length=0)

    fig.suptitle("Shift from wild type in each structural metric",
                 fontsize=10, y=1.0)

    out = args.output_dir
    (out / "plots").mkdir(parents=True, exist_ok=True)
    (out / "tables").mkdir(parents=True, exist_ok=True)
    f = out / "plots" / "figure3_structural_heatmap.png"
    fig.savefig(f, dpi=300)
    plt.close(fig)

    raw_holo.round(3).to_csv(out / "tables" / "figure3_holo_shift_vs_wt.csv")
    raw_apo.round(3).to_csv(out / "tables" / "figure3_apo_shift_vs_wt.csv")

    z_holo, z_apo = raw_holo, raw_apo
    print(f"holo metrics: {z_holo.shape[1]}, apo metrics: {z_apo.shape[1]}, "
          f"genotypes: {len(order)}")
    print("largest |shift| per holo metric, in its own units:")
    for c in raw_holo.columns:
        print(f"   {c:20s} {np.nanmax(np.abs(raw_holo[c].to_numpy())):8.2f} "
              f"{UNITS.get(c, '')}")
    print(f"\nwrote {f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
