#!/usr/bin/env python3
"""Does the NNIBP still collapse when the cavity is hydrated?

Two matched sets of short, densely sampled apo runs: one from the systems as
prepared, with the cavity left empty when DOR was deleted, and one with the
cavity filled with water, settled under protein restraints, then released.
The empty set plateaus near 70 A^3, about 30% below the unliganded crystal
structures (1DLO 101.7, 1HMV 103.8), which is the signature of relaxing into
a vacuum. This compares the two on the same axis.

Usage
-----
    PYTHONPATH=src python -m nnrti.cli.compare_collapse_hydration \
        --empty results/analysis/apo_collapse/apo_collapse_timecourse.csv \
        --filled results/analysis/apo_collapse/apo_collapse_solvated_timecourse.csv
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
APO_CRYSTAL = (101.7, 103.8)


def summarise(d: pd.DataFrame, label: str) -> pd.DataFrame:
    rows = []
    for (g, rep), s in d.groupby(["genotype", "replicate"]):
        s = s.sort_values("time_ps")
        v, t = s.pocket_volume.to_numpy(), s.time_ps.to_numpy()
        below = np.where(v <= max(APO_CRYSTAL))[0]
        rows.append({"set": label, "genotype": g, "replicate": rep,
                     "v0": v[0], "plateau": float(np.median(v[-100:])),
                     "t_cross_ps": float(t[below[0]]) if below.size else np.nan})
    return pd.DataFrame(rows)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--empty", type=Path,
                    default=ROOT / "results/analysis/apo_collapse/apo_collapse_timecourse.csv")
    ap.add_argument("--filled", type=Path,
                    default=ROOT / "results/analysis/apo_collapse/"
                                   "apo_collapse_solvated_timecourse.csv")
    ap.add_argument("--output-dir", type=Path,
                    default=ROOT / "results/analysis/apo_collapse")
    args = ap.parse_args(argv)

    out = []
    for p, label in ((args.empty, "empty cavity"), (args.filled, "water filled")):
        if not p.exists():
            print(f"missing {p}")
            continue
        out.append(summarise(pd.read_csv(p), label))
    if len(out) < 2:
        raise SystemExit("need both sets")
    s = pd.concat(out)

    print(f"{'set':14s} {'n':>3s} {'start':>7s} {'plateau':>9s} {'t to <=104':>11s}")
    for label, sub in s.groupby("set", sort=False):
        print(f"{label:14s} {len(sub):3d} {sub.v0.mean():7.1f} "
              f"{sub.plateau.mean():9.1f} {sub.t_cross_ps.median():8.0f} ps")
    print(f"\nunliganded crystals: {APO_CRYSTAL[0]} (1DLO), {APO_CRYSTAL[1]} (1HMV)")

    e = s[s.set == "empty cavity"].set_index(["genotype", "replicate"]).plateau
    f = s[s.set == "water filled"].set_index(["genotype", "replicate"]).plateau
    k = sorted(set(e.index) & set(f.index))
    if k:
        diff = f[k] - e[k]
        print(f"\nmatched replicates: {len(k)}")
        print(f"plateau difference (filled - empty): mean {diff.mean():+.1f} A^3, "
              f"median {diff.median():+.1f}, range {diff.min():+.0f} to {diff.max():+.0f}")
        print(f"filled plateau within the crystal band (90-115 A^3): "
              f"{int(((f[k] >= 90) & (f[k] <= 115)).sum())}/{len(k)}")
        print(f"empty  plateau within the crystal band: "
              f"{int(((e[k] >= 90) & (e[k] <= 115)).sum())}/{len(k)}")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    s.to_csv(args.output_dir / "collapse_hydration_comparison.csv", index=False)
    print(f"\nwrote {args.output_dir.relative_to(ROOT)}/collapse_hydration_comparison.csv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
