#!/usr/bin/env python3
"""Protein-level changes: per-residue RMSF and secondary structure, WT-referenced.

Reviewer 1 (comment 4) asked for protein-level observables that the manuscript
did not show -- RMSF, secondary-structure changes and pocket conformational
distributions -- to separate direct changes in DOR contacts from broader
remodelling of the pocket.

For each holo replicate this computes:

  * per-residue Ca RMSF after superposition on the pocket-lining residues, so
    the fluctuation reported is motion relative to the pocket rather than
    whole-molecule tumbling;
  * per-residue secondary structure (DSSP, three-state) as the fraction of
    frames in helix and in strand.

Both are then referenced to WT. The regions reported in the summary are the
ones the reviewers raised: the NNIBP lining, the 225-235 primer grip, and the
beta6-beta10-beta9 sheet around 100-110 and 179-190.

Usage
-----
    PYTHONPATH=src python -m nnrti.cli.compute_protein_level_changes
"""

from __future__ import annotations

import argparse
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

RESID_OFFSET = -3
NNIBP_P66 = (100, 101, 103, 106, 107, 108, 179, 181, 188, 189, 190, 227, 229, 234, 318)
REGIONS = {
    "NNIBP lining": NNIBP_P66,
    "primer grip 225-235": tuple(range(225, 236)),
    "beta6 100-110": tuple(range(100, 111)),
    "beta9-beta10 179-190": tuple(range(179, 191)),
}


def main(argv: list[str] | None = None) -> int:
    import mdtraj as md

    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--runs-root", type=Path, default=Path("results/md_runs"))
    p.add_argument("--output-dir", type=Path, default=Path("results/analysis/protein_level"))
    p.add_argument("--max-frames", type=int, default=200,
                   help="Subsample to at most this many frames per replicate.")
    args = p.parse_args(argv)

    dirs = sorted(d for d in args.runs_root.glob("*/rep_*")
                  if "apo" not in d.parts and (list(d.glob("*analysis.dcd"))))
    if not dirs:
        raise SystemExit("no holo replicates found")

    rmsf_rows, ss_rows = [], []
    for d in dirs:
        label, rep = d.parent.name, int(d.name.split("_")[-1])
        dcd = sorted(d.glob("*analysis.dcd"))[0]
        top = sorted(d.glob("*analysis_topology.pdb"))[0]
        try:
            t = md.load(str(dcd), top=str(top))
        except Exception as e:
            print(f"  {label} rep{rep}: load failed ({e})", flush=True)
            continue
        if t.n_frames > args.max_frames:
            t = t[np.linspace(0, t.n_frames - 1, args.max_frames).astype(int)]

        pocket_ca = [a.index for a in t.topology.atoms
                     if a.name == "CA" and a.residue.resSeq in
                     [r + RESID_OFFSET for r in NNIBP_P66]]
        if not pocket_ca:
            print(f"  {label} rep{rep}: no pocket CA found", flush=True)
            continue
        t.superpose(t, 0, atom_indices=pocket_ca)

        ca = [a for a in t.topology.atoms if a.name == "CA"]
        idx = [a.index for a in ca]
        xyz = t.xyz[:, idx, :] * 10.0
        rmsf = np.sqrt(((xyz - xyz.mean(0)) ** 2).sum(-1).mean(0))
        for a, v in zip(ca, rmsf):
            rmsf_rows.append({"safe_label": label, "replicate": rep,
                              "resSeq": a.residue.resSeq,
                              "canonical": a.residue.resSeq - RESID_OFFSET,
                              "chain": a.residue.chain.index, "rmsf": float(v)})

        try:
            dssp = md.compute_dssp(t, simplified=True)
            residues = [r for r in t.topology.residues if r.is_protein]
            for j, r in enumerate(residues):
                col = dssp[:, j]
                ss_rows.append({"safe_label": label, "replicate": rep,
                                "resSeq": r.resSeq, "canonical": r.resSeq - RESID_OFFSET,
                                "chain": r.chain.index,
                                "frac_helix": float((col == "H").mean()),
                                "frac_strand": float((col == "E").mean())})
        except Exception as e:
            print(f"  {label} rep{rep}: DSSP failed ({e})", flush=True)
        print(f"  {label:14s} rep{rep}  {t.n_frames} frames  "
              f"mean pocket RMSF {rmsf[[ca.index(a) for a in ca if a.index in pocket_ca]].mean():.2f} A",
              flush=True)

    out = args.output_dir / "tables"
    out.mkdir(parents=True, exist_ok=True)
    rmsf = pd.DataFrame(rmsf_rows)
    rmsf.to_csv(out / "ca_rmsf_per_residue.csv", index=False)
    ss = pd.DataFrame(ss_rows)
    if not ss.empty:
        ss.to_csv(out / "secondary_structure_per_residue.csv", index=False)

    # WT-referenced regional summary
    gm = rmsf.groupby(["safe_label", "chain", "canonical"]).rmsf.mean().reset_index()
    wt = gm[gm.safe_label == "wt"].set_index(["chain", "canonical"]).rmsf
    rows = []
    for lbl, g in gm[gm.safe_label != "wt"].groupby("safe_label"):
        g = g.set_index(["chain", "canonical"])
        delta = (g.rmsf - wt).dropna()
        row = {"safe_label": lbl}
        for name, resids in REGIONS.items():
            sel = delta[[i for i in delta.index if i[1] in resids and i[0] == 0]]
            row[name] = round(float(sel.mean()), 3) if len(sel) else float("nan")
        rows.append(row)
    summ = pd.DataFrame(rows).set_index("safe_label").sort_index()
    summ.to_csv(out / "rmsf_delta_by_region.csv")

    print("\n=== mean Ca RMSF change vs WT, by region (A) ===")
    print(summ.to_string())
    print(f"\nwrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
