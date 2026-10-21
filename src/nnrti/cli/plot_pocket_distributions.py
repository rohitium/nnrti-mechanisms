#!/usr/bin/env python3
"""NNIBP volume distributions, bound and unliganded.

Reviewer 1 (comment 4) asked for pocket conformational distributions alongside
RMSF and secondary structure, to separate changes in DOR contacts from broader
conformational effects on the pocket.

Means hide the thing that matters here. A pocket that is bimodal, or one whose
spread doubles without its mean moving, reads identically to a stable pocket in
Table 3. Plotting the per-frame distribution shows which it is.

Both states are shown on one axis because the comparison that matters is
between them: the NNIBP is an induced pocket, so the gap between the bound and
unliganded distributions is the conformational change DOR has to pay for, and
whether a mutation narrows or widens that gap is the mechanistic question
behind Reviewer 1's argument.

The unliganded panel uses the matched 10 ns window, for the reason given in
compare_apo_ensembles: WT apo is 10 ns and is the reference.

Usage
-----
    PYTHONPATH=src python -m nnrti.cli.plot_pocket_distributions
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

CATS = {"Susceptible": "#2c7fb8", "Resistant": "#d7301f", "Uncertain": "#756bb1"}


def main(argv: list[str] | None = None) -> int:
    from nnrti.analysis.panel import RESISTANT, SUSCEPTIBLE, UNCERTAIN

    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--holo", type=Path,
                   default=Path("results/.checkpoints/.checkpoint_pocket_volume_profiles.csv"))
    p.add_argument("--apo", type=Path,
                   default=Path("results/analysis/apo_ensembles/tables/apo_frame_metrics.csv"))
    p.add_argument("--output-dir", type=Path, default=Path("results/analysis/mechanisms/plots"))
    p.add_argument("--dpi", type=int, default=300)
    args = p.parse_args(argv)

    holo = pd.read_csv(args.holo)
    holo = holo.rename(columns={"pocket_volume_proxy_angstrom3": "volume"})
    apo = pd.read_csv(args.apo)
    apo = apo[apo.window == "matched_10ns"].rename(columns={"pocket_volume": "volume"})
    # The holo tables capitalise the genotype label and the apo runs do not;
    # normalise before any join or the two states never line up.
    holo["safe_label"] = holo["safe_label"].str.lower()
    apo["safe_label"] = apo["safe_label"].str.lower()

    def cat(label: str) -> str:
        g = label.replace("_", "+")
        if g.lower() == "wt":
            return "Wild type"
        for grp, name in ((SUSCEPTIBLE, "Susceptible"), (RESISTANT, "Resistant"),
                          (UNCERTAIN, "Uncertain")):
            if g.upper() in {x.upper() for x in grp}:
                return name
        return "Uncertain"

    order = sorted(set(holo.safe_label) | set(apo.safe_label),
                   key=lambda s: (s != "wt", cat(s), s))

    fig, axes = plt.subplots(1, 2, figsize=(13.0, 7.4), sharey=True)
    for ax, (d, title, state) in zip(axes, [
            (holo, "A.  Bound (holo), full production", "holo"),
            (apo, "B.  Unliganded (apo), matched 10 ns window", "apo")]):
        for i, lbl in enumerate(order):
            v = d.loc[d.safe_label == lbl, "volume"].dropna().to_numpy()
            if v.size == 0:
                continue
            c = "black" if lbl == "wt" else CATS.get(cat(lbl), "0.4")
            parts = ax.violinplot([v], positions=[i], vert=False, widths=0.85,
                                  showextrema=False, showmedians=True)
            for b in parts["bodies"]:
                b.set_facecolor(c); b.set_alpha(0.55); b.set_edgecolor(c); b.set_linewidth(0.6)
            parts["cmedians"].set_color("black"); parts["cmedians"].set_linewidth(1.1)
        # Labels come from `order`, not from what this panel happens to contain:
        # the axes share a y-axis, so a per-panel list would be overwritten by
        # whichever panel was drawn last and silently mislabel the other.
        ax.set_yticks(range(len(order)))
        ax.set_yticklabels([l.replace("_", "+") for l in order], fontsize=8)
        if "wt" in order:
            wtv = d.loc[d.safe_label == "wt", "volume"].median()
            if np.isfinite(wtv):
                ax.axvline(wtv, color="black", ls="--", lw=0.9, alpha=0.7)
                ax.text(wtv, len(order) - 0.3, " WT median", fontsize=7, color="0.3")
        # A handful of frames carry non-physical volumes (the proxy losing its
        # footing, as in the G190E apo runs); they would otherwise set the scale
        # and flatten every distribution into a line.
        hi = float(np.nanpercentile(d["volume"].to_numpy(), 99.5))
        n_clipped = int((d["volume"] > hi).sum())
        ax.set_xlim(0, hi * 1.08)
        if n_clipped:
            ax.text(0.99, 0.01, f"{n_clipped} frames beyond axis", transform=ax.transAxes,
                    fontsize=6.5, color="0.45", ha="right", va="bottom")
        ax.set_xlabel("NNIBP volume (Å$^3$)", fontsize=9)
        ax.set_title(title, fontsize=10, loc="left")
        ax.grid(axis="x", alpha=0.25, lw=0.5)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
    axes[0].invert_yaxis()

    handles = [plt.Line2D([], [], marker="s", ls="", color=c, markersize=8, label=k)
               for k, c in CATS.items()]
    handles.insert(0, plt.Line2D([], [], marker="s", ls="", color="black", markersize=8,
                                 label="Wild type"))
    axes[1].legend(handles=handles, fontsize=7.5, loc="lower right", frameon=False)
    fig.tight_layout()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(args.output_dir / f"pocket_volume_distributions.{ext}",
                    dpi=args.dpi, bbox_inches="tight")

    hs = holo.groupby("safe_label").volume.agg(["median", "std"]).round(1)
    as_ = apo.groupby("safe_label").volume.agg(["median", "std"]).round(1)
    j = hs.join(as_, lsuffix="_holo", rsuffix="_apo", how="outer")
    j["induced_gap"] = (j["median_holo"] - j["median_apo"]).round(1)
    j.to_csv(args.output_dir.parent / "tables" / "pocket_volume_distributions.csv")
    print(f"wrote {args.output_dir}/pocket_volume_distributions.png|pdf\n")
    print(j.dropna(subset=["median_apo"]).to_string())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
