#!/usr/bin/env python3
"""Cost of opening the NNIBP: the bound pocket volume minus the unliganded one.

Reviewer 1 (comment 3) asks whether mutation-induced changes occur in both apo
and holo RT, since a change present in both largely cancels in a thermodynamic
cycle. The opening cost is the simplest expression of that: how much the pocket
has to grow for DOR to bind, per genotype, and how far that differs from wild
type.

Both terms are means over replicates, so the difference carries both errors.
Earlier versions of this table reported the cost without one, which made an
8 A^3 difference and a 99 A^3 difference look like the same kind of statement.

Usage
-----
    PYTHONPATH=src python -m nnrti.cli.compute_pocket_opening_cost
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def _display(label: str) -> str:
    return "WT" if label.lower() == "wt" else label.upper().replace("_", "+")


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--apo-frames", type=Path,
                   default=Path("results/analysis/apo_ensembles/tables/"
                                "apo_frame_metrics.csv"))
    p.add_argument("--holo-per-rep", type=Path,
                   default=Path("results/analysis/pocket_volume/tables/"
                                "pocket_volume_per_rep.csv"))
    p.add_argument("--pocket-rg-max", type=float, default=10.0)
    p.add_argument("--output", type=Path,
                   default=Path("results/analysis/apo_ensembles/tables/"
                                "pocket_opening_cost.csv"))
    args = p.parse_args(argv)

    # Apo: replicate means over the matched window, excluding frames whose
    # pocket has come apart and replicates too short for the window.
    f = pd.read_csv(args.apo_frames)
    # Same three filters compare_apo_ensembles applies, in the same order:
    # intact pocket, the matched window, and replicates long enough to have
    # one. Leaving the last out put apo v106a rep_02 (34 ns) into a 100 ns
    # comparison and moved V106A's apo volume by 24 A^3.
    if "pocket_intact" in f.columns:
        f = f[f.pocket_intact]
    else:
        f = f[f.pocket_rg <= args.pocket_rg_max]
    f = f[f.window == "matched"]
    total = f.groupby(["safe_label", "replicate"]).total_ns.first()
    window_ns = f.time_ns.max()
    keep = {k for k, v in total.items() if v >= window_ns * 0.98}
    dropped = sorted(set(total.index) - keep)
    if dropped:
        print("excluded, shorter than the matched window: "
              + ", ".join(f"{l} rep{r} ({total[(l, r)]:.1f} ns)" for l, r in dropped))
    f = f[[(l, r) in keep for l, r in zip(f.safe_label, f.replicate)]]
    apo_rep = f.groupby(["safe_label", "replicate"], as_index=False).pocket_volume.mean()
    apo = apo_rep.groupby("safe_label").pocket_volume.agg(["mean", "sem", "size"])
    apo.columns = ["apo_volume_A3", "apo_sem_A3", "n_apo"]
    apo.index = [_display(i) for i in apo.index]

    holo_rep = pd.read_csv(args.holo_per_rep)
    holo = holo_rep.groupby("mutation").pocket_volume_mean.agg(["mean", "sem", "size"])
    holo.columns = ["holo_volume_A3", "holo_sem_A3", "n_holo"]

    t = holo.join(apo, how="inner")
    t["opening_A3"] = t.holo_volume_A3 - t.apo_volume_A3
    t["opening_sem_A3"] = np.sqrt(t.holo_sem_A3 ** 2 + t.apo_sem_A3 ** 2)

    wt = t.loc["WT"]
    t["vs_WT_A3"] = t.opening_A3 - wt.opening_A3
    # Against wild type, four means are subtracted, so four errors add.
    t["vs_WT_sem_A3"] = np.sqrt(t.opening_sem_A3 ** 2 + wt.opening_sem_A3 ** 2)
    t["vs_WT_sem_ratio"] = (t.vs_WT_A3.abs() / t.vs_WT_sem_A3).replace(np.inf, np.nan)

    t = t.round(2).sort_values("opening_A3", ascending=False)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    t.to_csv(args.output)

    print(t[["n_apo", "apo_volume_A3", "apo_sem_A3", "holo_volume_A3",
             "opening_A3", "opening_sem_A3", "vs_WT_A3", "vs_WT_sem_A3",
             "vs_WT_sem_ratio"]].to_string())
    print(f"\nWT opening cost: {wt.opening_A3:.1f} +/- {wt.opening_sem_A3:.1f} A^3")
    res = t[(t.index != "WT") & (t.vs_WT_sem_ratio > 2.306)]
    print(f"genotypes whose opening cost differs from WT at 95% "
          f"(>2.31 SEM, df~8): {', '.join(res.index) if len(res) else 'none'}")
    print(f"\nwrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
