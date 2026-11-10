#!/usr/bin/env python3
"""Per-frame NNIBP volume for the bound state, on the apo convention.

The apo per-frame volumes already exist, in
results/analysis/apo_ensembles/tables/apo_frame_metrics.csv, and Supplementary
Table 4's unliganded column is their mean. The bound column comes from
compute_pocket_volume, which writes a mean and a standard deviation but not
the frames, so there is nothing to draw a bound distribution over.

That script also cannot simply be re-run: it reads the *_analysis_pbcfix.dcd
sidecar, and no holo replicate has one on this machine -- DCDs are gitignored,
so the sidecars the published numbers were made from are not here. It would
skip all 60 replicates.

So the correction is done in memory instead, with the convention
compare_apo_ensembles uses for the apo volumes: make_molecules_whole, then
image_molecules anchored on all protein atoms, then the coordinates injected
frame by frame into a topology-only Universe and scored by the same
pocket_volume_proxy_from_universe. compute_protein_level_changes already
establishes that in-memory imaging reproduces a pbcfix sidecar exactly, and
--validate below is the direct check: the mean over the same 40-frame subset
compute_pocket_volume scored must come back as the published value.

ligand_resname is ignored by the proxy, so this is the protein's own cavity
with DOR not counted as an occluder -- the same quantity in both states, which
is what makes the bound/unliganded comparison meaningful at all.

Usage
-----
    PYTHONPATH=src python -m nnrti.cli.compute_holo_volume_frames
    PYTHONPATH=src python -m nnrti.cli.compute_holo_volume_frames --validate
"""
from __future__ import annotations

import argparse
import logging
from pathlib import Path

import numpy as np
import pandas as pd

from nnrti.analysis.metrics import pocket_volume_proxy_from_universe
from nnrti.cli.compute_block_metrics import _display, holo_frame_times_ns
from nnrti.cli.compute_pocket_volume import _frame_subset
from nnrti.cli.compute_protein_level_changes import NNIBP_P66, RESID_OFFSET
from nnrti.cli.find_pocket_surrogates import REPO_ROOT

#: The frame count compute_pocket_volume scored per replicate, whose mean is
#: the published bound volume. Used by --validate only.
PUBLISHED_MAX_FRAMES = 40


def main(argv: list[str] | None = None) -> int:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    import MDAnalysis as mda
    import mdtraj as md

    from nnrti.analysis.metrics import _p66_segid

    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--max-frames", type=int, default=200,
                   help="Frames scored per replicate for the distribution.")
    p.add_argument("--genotypes", nargs="*", default=None)
    p.add_argument("--validate", action="store_true",
                   help="Also score the published 40-frame subset and compare "
                        "its mean against pocket_volume_per_rep.csv.")
    p.add_argument("--output", type=Path,
                   default=REPO_ROOT / "results/analysis/pocket_volume/tables/"
                                       "holo_volume_per_frame.csv")
    args = p.parse_args(argv)

    root = Path("results/md_runs")
    dirs = sorted(d for d in root.glob("*/rep_*")
                  if "apo" not in d.parts and "pose" not in str(d)
                  and "solvent" not in str(d) and "collapse" not in str(d)
                  and list(d.glob("*analysis.dcd")))
    if args.genotypes:
        want = {g.lower() for g in args.genotypes}
        dirs = [d for d in dirs if d.parent.name.lower() in want]

    rows, checks = [], []
    for d in dirs:
        label, rep = d.parent.name, int(d.name.split("_")[-1])
        try:
            dcd = sorted(q for q in d.glob("*analysis.dcd") if "pbcfix" not in q.name)[0]
            top = sorted(q for q in d.glob("*analysis_topology.pdb")
                         if "from_md" not in q.name)[0]
            t = md.load(str(dcd), top=str(top))
            t.make_molecules_whole(inplace=True)
            prot = t.topology.select("protein")
            if len(prot):
                t.image_molecules(inplace=True,
                                  anchor_molecules=[{t.topology.atom(i) for i in prot}])
            times = holo_frame_times_ns(d, t.n_frames)
            xyz = t.xyz * 10.0

            u = mda.Universe(str(top))
            seg = _p66_segid(u)
            resid_sel = " ".join(str(r + RESID_OFFSET) for r in NNIBP_P66)
            ca = u.select_atoms(f"segid {seg} and name CA and resid {resid_sel}")
            if ca.n_atoms != len(NNIBP_P66):
                raise ValueError(f"{ca.n_atoms} pocket CA atoms, "
                                 f"expected {len(NNIBP_P66)}")

            def score(idx: np.ndarray) -> np.ndarray:
                out = np.empty(len(idx))
                for k, fi in enumerate(idx):
                    u.atoms.positions = xyz[int(fi)]
                    out[k] = float(pocket_volume_proxy_from_universe(u, ""))
                return out

            n = t.n_frames
            idx = _frame_subset(n, args.max_frames, max(1, n // args.max_frames))
            vals = score(idx)
            for fi, v in zip(idx, vals):
                rows.append({"mutation": _display(label), "replicate": rep,
                             "frame": int(fi), "time_ns": round(float(times[fi]), 3),
                             "pocket_volume": float(v),
                             "total_ns": round(float(times[-1]), 2)})

            if args.validate:
                pidx = _frame_subset(n, PUBLISHED_MAX_FRAMES,
                                     max(1, n // PUBLISHED_MAX_FRAMES))
                checks.append({"mutation": _display(label), "replicate": rep,
                               "n_frames": len(pidx),
                               "recomputed_mean": float(score(pidx).mean())})
            logging.info("%-12s rep%02d  %d of %d frames  mean %.1f A^3",
                         label, rep, len(idx), n, vals.mean())
        except Exception as exc:
            logging.error("FAILED %s rep%02d: %s", label, rep, exc)

    if not rows:
        return 1
    args.output.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(args.output, index=False)
    logging.info("wrote %s (%d frames)", args.output, len(rows))

    if checks:
        pub = pd.read_csv(REPO_ROOT / "results/analysis/pocket_volume/tables/"
                                      "pocket_volume_per_rep.csv")
        cmp = pd.DataFrame(checks).merge(
            pub[["mutation", "replicate", "pocket_volume_mean", "pocket_volume_std"]],
            on=["mutation", "replicate"], how="left")
        cmp["diff"] = cmp.recomputed_mean - cmp.pocket_volume_mean
        cmp["diff_in_sd"] = cmp["diff"] / cmp.pocket_volume_std
        print(cmp.round(3).to_string(index=False))
        print(f"\nmax |diff| {cmp['diff'].abs().max():.2f} A^3, "
              f"max |diff| in published SD {cmp['diff_in_sd'].abs().max():.2f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
