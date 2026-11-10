#!/usr/bin/env python3
"""Per-frame chi1 of Tyr181 and Tyr188, in the apo and the holo state.

compare_chi1_apo_holo.py reports one circular mean per replicate, which is all
the heatmap in Figure 5 needs. A violin needs the frames themselves: three
replicate means cannot carry a density, and a rotamer that is bimodal reads
identically to a rigid side chain once it has been averaged.

Every frame of the production trajectory is emitted, which is the window
compare_chi1_apo_holo averages over and therefore the window behind
Supplementary Table 4's chi1 columns. The circular mean of what this writes
reproduces that table exactly, so the violin and the published number are the
same quantity.

The accompanying time axis is state-dependent, via
compute_block_metrics.frame_times. It must not come from
compare_apo_ensembles.frame_times_ns alone: that reconstructs the apo
two-stride schedule, and 20 of the 60 holo replicates hold more than 200
frames, for which it returns anything from 29 to 180 ns against a 100 ns step
record. The axis is not used to select frames here, only to label them.

A torsion is intra-residue, so no periodic-image correction is needed in
principle; atoms of one residue are nonetheless imaged against the residue's
own N atom, because a residue straddling a box face would otherwise give a
geometrically meaningless angle with nothing to flag it.

Genotypes that remove the residue are excluded from its series: Y181C has no
Tyr181, Y188L no Tyr188.

Usage
-----
    PYTHONPATH=src python -m nnrti.cli.compute_chi1_frame_series
"""
from __future__ import annotations

import argparse
import logging
from pathlib import Path

import numpy as np
import pandas as pd

from nnrti.cli.compare_chi1_apo_holo import RESID_OFFSET, RESIDUES, circular_mean
from nnrti.cli.compute_block_metrics import frame_times
from nnrti.cli.compute_mechanism_coordinates import _replicate_inputs
from nnrti.cli.find_pocket_surrogates import (
    REPO_ROOT,
    _apo_replicates,
    _mic,
    _run_ns,
)

CHI1 = ("N", "CA", "CB", "CG")


def _dihedrals(p: np.ndarray) -> np.ndarray:
    """Signed dihedral in degrees for a stack of four-point sets, (n, 4, 3)."""
    b0 = p[:, 0] - p[:, 1]
    b1 = p[:, 2] - p[:, 1]
    b2 = p[:, 3] - p[:, 2]
    b1 = b1 / np.linalg.norm(b1, axis=1, keepdims=True)
    v = b0 - (b0 * b1).sum(1, keepdims=True) * b1
    w = b2 - (b2 * b1).sum(1, keepdims=True) * b1
    return np.degrees(np.arctan2((np.cross(b1, v) * w).sum(1), (v * w).sum(1)))


def _chi1_atoms(top, resid: int) -> list[int] | None:
    """The four chi1 atom indices of p66 residue `resid`, if it is a tyrosine.

    Leucine also has N-CA-CB-CG, so an atom-presence test alone would let
    Y188L through and compare a leucine chi1 against a tyrosine's.
    """
    res = [r for r in top.residues
           if r.resSeq == resid + RESID_OFFSET and r.chain.index == 0]
    if len(res) != 1 or res[0].name != "TYR":
        return None
    by_name = {a.name: a.index for a in res[0].atoms}
    if not set(CHI1) <= set(by_name):
        return None
    return [by_name[n] for n in CHI1]


def main(argv: list[str] | None = None) -> int:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    import mdtraj as md

    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--manifest", type=Path, default=Path("manifests/md_manifest.csv"))
    p.add_argument("--min-ns", type=float, default=80.0)
    p.add_argument("--mutations", nargs="*", default=None)
    p.add_argument("--output", type=Path,
                   default=REPO_ROOT / "results/analysis/chi1_apo_holo/tables/"
                                       "chi1_per_frame.csv")
    args = p.parse_args(argv)

    mf = pd.read_csv(args.manifest)
    if args.mutations:
        mf = mf[mf["mutation"].isin(args.mutations)]

    jobs = []
    for _, r in mf.iterrows():
        try:
            topo, dcd, _ = _replicate_inputs(r, REPO_ROOT)
            jobs.append(("holo", str(r["mutation"]), int(r["replicate"]), topo, dcd))
        except Exception as exc:
            logging.error("holo %s rep%s: %s", r["mutation"], r["replicate"], exc)
    for mutation in dict.fromkeys(str(m) for m in mf["mutation"]):
        for rep, topo, dcd in _apo_replicates(mutation):
            jobs.append(("apo", mutation, rep, topo, dcd))

    rows = []
    for state, mutation, rep, topo, dcd in jobs:
        try:
            total_ns = _run_ns(Path(dcd).parent)
            if total_ns < args.min_ns:
                logging.warning("skipped %s %s rep%02d: only %.0f ns on disk",
                                mutation, state, rep, total_ns)
                continue
            traj = md.load(str(dcd), top=str(topo))
            box = (traj.unitcell_lengths * 10.0) if traj.unitcell_lengths is not None else None
            xyz = traj.xyz * 10.0
            times = frame_times(state, Path(dcd).parent, traj.n_frames)
            for label, canonical in RESIDUES.items():
                idx = _chi1_atoms(traj.topology, canonical)
                if idx is None:
                    continue
                pts = xyz[:, idx, :]
                if box is not None:
                    ref = pts[:, 0:1, :]
                    for k in range(pts.shape[0]):
                        pts[k] = ref[k] + _mic(pts[k] - ref[k], box[k])
                ang = _dihedrals(pts)
                full_mean = circular_mean(ang)
                for f, (a, tt) in enumerate(zip(ang, times)):
                    rows.append(dict(mutation=mutation, state=state, replicate=rep,
                                     residue=label, frame=int(f),
                                     time_ns=round(float(tt), 3),
                                     run_ns=round(total_ns, 1),
                                     chi1_deg=round(float(a), 3),
                                     chi1_full_traj_mean=round(full_mean, 2)))
            logging.info("%s %s rep%02d: %d frames over %.1f ns",
                         mutation, state, rep, traj.n_frames, times[-1])
        except Exception as exc:
            logging.error("FAILED %s %s rep%s: %s", mutation, state, rep, exc)

    if not rows:
        return 1
    df = pd.DataFrame(rows)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.output, index=False)
    logging.info("wrote %s (%d rows)", args.output, len(df))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
