#!/usr/bin/env python3
"""Pocket volume against time for the densely sampled apo runs.

The runs save protein coordinates every picosecond. Volumes are computed here,
locally, with the same routine used for every other pocket number in the paper,
so the collapse curves and the ensemble values in Table 7 are on one scale.

Usage
-----
    PYTHONPATH=src python -m nnrti.cli.analyse_apo_collapse --coords DIR
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]


def main(argv: list[str] | None = None) -> int:
    import MDAnalysis as mda
    from nnrti.analysis.metrics import pocket_volume_proxy_from_universe

    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--coords", type=Path, required=True)
    ap.add_argument("--apo-root", type=Path, default=ROOT / "results/md_runs/apo")
    ap.add_argument("--frame-ps", type=float, default=1.0)
    ap.add_argument("--output-dir", type=Path,
                    default=ROOT / "results/analysis/apo_collapse")
    args = ap.parse_args(argv)

    files = sorted(args.coords.glob("*_rep*_coords.npy"))
    if not files:
        raise SystemExit(f"no *_coords.npy under {args.coords}")

    rows = []
    for f in files:
        m = re.match(r"(.+)_rep(\d+)_coords\.npy$", f.name)
        if not m:
            continue
        g, rep = m.group(1), int(m.group(2))
        idx_f = f.with_name(f"{g}_rep{rep:02d}_atomindex.npy")
        start = args.apo_root / g / f"rep_{rep:02d}" / "assets" / \
            f"{g}_apo_md_rep{rep:02d}_start.pdb"
        if not (idx_f.exists() and start.exists()):
            print(f"  skip {g} rep{rep}: missing index or start structure")
            continue
        coords = np.load(f)
        idx = np.load(idx_f)
        u = mda.Universe(str(start))
        if coords.shape[1] != idx.size:
            print(f"  skip {g} rep{rep}: {coords.shape[1]} coords vs {idx.size} indices")
            continue
        for i in range(coords.shape[0]):
            u.atoms[idx].positions = coords[i]
            rows.append({"genotype": g, "replicate": rep,
                         "time_ps": round(i * args.frame_ps, 3),
                         "pocket_volume": round(
                             float(pocket_volume_proxy_from_universe(u, "")), 2)})
        v = [r["pocket_volume"] for r in rows if r["genotype"] == g
             and r["replicate"] == rep]
        print(f"  {g:14s} rep{rep}: t=0 {v[0]:6.1f}  t=5 {v[5]:6.1f}  "
              f"final {v[-1]:6.1f}  median(last 100) {np.median(v[-100:]):6.1f}")

    df = pd.DataFrame(rows)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    out = args.output_dir / "apo_collapse_timecourse.csv"
    df.to_csv(out, index=False)
    print(f"\nwrote {out}  ({len(df)} rows, "
          f"{df.genotype.nunique()} genotypes)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
