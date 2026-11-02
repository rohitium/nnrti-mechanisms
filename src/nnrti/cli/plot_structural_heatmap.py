#!/usr/bin/env python3
"""Figure 3: WT-referenced structural metrics across the panel.

Reviewer 1 (comment 4) asked for a WT-referenced heatmap summarising the
structural metrics, because they are currently spread across dense tables and
patterns in DOR contacts, pocket size and remodelling are hard to see.

Each cell is a genotype's shift from wild type in units of the pooled
replicate standard deviation of that metric. Raw units are not comparable --
an angstrom of contact distance and a cubic angstrom of pocket volume are not
the same size of change -- so expressing every metric as "how many replicate
error bars from WT" is what makes one colour scale legitimate across columns.
The per-genotype values in their own units are in Supplementary Table 4.

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

#: Colour saturates here, in replicate standard deviations.
CLIP = 4.0


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
                  metrics: dict[str, str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Shift from WT in pooled replicate SDs, plus the raw shift."""
    d = per_rep[per_rep.structure_type == state]
    z, raw = {}, {}
    for label, col in metrics.items():
        if col not in d or d[col].isna().all():
            continue
        g = d.groupby("genotype")[col]
        means = g.mean()
        # Pool the within-genotype replicate SD across the panel: a single
        # scale for the metric, rather than one noisy SD per genotype.
        sd = float(np.sqrt((g.var(ddof=1).dropna()).mean()))
        if not np.isfinite(sd) or sd == 0:
            continue
        wt = means.get("WT")
        if wt is None or not np.isfinite(wt):
            continue
        raw[label] = (means - wt).round(3)
        z[label] = ((means - wt) / sd).round(3)
    return pd.DataFrame(z), pd.DataFrame(raw)


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
    z_holo, raw_holo = wt_referenced(per_rep, "holo", all_metrics)
    apo_metrics = {**POCKET_METRICS, **PROTEIN_METRICS}
    z_apo, raw_apo = wt_referenced(per_rep, "apo", apo_metrics)

    order = sorted((g for g in z_holo.index if g != "WT"),
                   key=lambda g: (np.nan_to_num(fold.get(g, np.nan), nan=-1)))
    z_holo = z_holo.loc[order]
    z_apo = z_apo.reindex(order)

    norm = TwoSlopeNorm(vmin=-CLIP, vcenter=0.0, vmax=CLIP)
    panels = [("DOR-bound (holo)", z_holo), ("Unliganded (apo)", z_apo)]
    widths = [p[1].shape[1] for p in panels]
    fig, axes = plt.subplots(
        1, 2, figsize=(0.30 * sum(widths) + 4.2, 0.26 * len(order) + 1.9),
        gridspec_kw={"width_ratios": widths, "wspace": 0.06})

    for ax, (title, z) in zip(axes, panels):
        im = ax.imshow(z.to_numpy(dtype=float), cmap="RdBu_r", norm=norm,
                       aspect="auto")
        ax.set_xticks(range(z.shape[1]))
        ax.set_xticklabels(z.columns, rotation=90, ha="center")
        ax.set_yticks(range(len(order)))
        ax.set_yticklabels([f"{g}  ({fold.get(g):g}×)" if np.isfinite(
            fold.get(g, np.nan)) else g for g in order] if ax is axes[0] else [])
        ax.set_title(title, fontsize=9)
        ax.set_xticks(np.arange(-0.5, z.shape[1], 1), minor=True)
        ax.set_yticks(np.arange(-0.5, len(order), 1), minor=True)
        ax.grid(which="minor", color="white", lw=0.6)
        ax.tick_params(which="minor", length=0)
        # A blank cell is a metric that cannot exist, not a missing measurement.
        for i in range(z.shape[0]):
            for j in range(z.shape[1]):
                if not np.isfinite(z.iat[i, j]):
                    ax.text(j, i, "·", ha="center", va="center",
                            color="0.6", fontsize=7)

    cb = fig.colorbar(im, ax=axes, fraction=0.025, pad=0.02, extend="both")
    cb.set_label("shift from WT (replicate SD)", fontsize=7)
    fig.suptitle("Structural change across the panel, referenced to wild type",
                 fontsize=10, y=1.0)

    out = args.output_dir
    (out / "plots").mkdir(parents=True, exist_ok=True)
    (out / "tables").mkdir(parents=True, exist_ok=True)
    f = out / "plots" / "figure3_structural_heatmap.png"
    fig.savefig(f, dpi=300)
    plt.close(fig)

    for name, df in (("holo_z", z_holo), ("apo_z", z_apo),
                     ("holo_raw", raw_holo.loc[order]),
                     ("apo_raw", raw_apo.reindex(order))):
        df.to_csv(out / "tables" / f"figure3_{name}.csv")

    print(f"holo metrics: {z_holo.shape[1]}, apo metrics: {z_apo.shape[1]}, "
          f"genotypes: {len(order)}")
    print(f"largest |shift| holo: {np.nanmax(np.abs(z_holo.to_numpy())):.1f} SD")
    print(f"largest |shift| apo : {np.nanmax(np.abs(z_apo.to_numpy())):.1f} SD")
    print(f"\nwrote {f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
