#!/usr/bin/env python3
"""NNIBP pocket volume per frame, in the apo and the holo state.

V(NNIBP) is the one Table 3 column that needs no surrogate: it is a property
of the protein alone, so the same quantity exists with and without doravirine
and the double difference

    ddV = [V(mut, holo) - V(WT, holo)] - [V(mut, apo) - V(WT, apo)]

can be formed directly. The published per-replicate volumes cover the holo
state only, and the apo volumes on disk predate the current apo runs, so both
states are recomputed here on one matched trailing window with one periodicity
convention.

Periodicity is handled by minimum image against the pocket centroid rather
than by MDTraj's make_molecules_whole. Bond-graph unwrapping propagates a
single over-long bond through the rest of the chain, which is what inflated
the apo G190E pocket past 2000 A^3 in the earlier pass. Only atoms within
about 12 A of the pocket centre can touch the grid, and nothing in RT sits
more than ~70 A from it against a 128 A box, so the image is unambiguous.

Usage
-----
    PYTHONPATH=src python -m nnrti.cli.compute_nnibp_volume_both_states
"""
from __future__ import annotations

import argparse
import logging
from pathlib import Path

import numpy as np
import pandas as pd

from nnrti.cli.compute_mechanism_coordinates import OFFSET, _replicate_inputs
from nnrti.cli.find_pocket_surrogates import (
    REPO_ROOT,
    _apo_replicates,
    _mic,
    _run_ns,
    _window_mask,
)

#: Cilento, Kirby & Sarafianos: the 15 p66 residues lining the NNIBP, plus
#: Glu138 of p51. Canonical numbering; the topologies run three lower.
NNIBP_P66 = (100, 101, 103, 106, 107, 108, 179, 181, 188, 189, 190, 227, 229, 234, 318)


def main(argv: list[str] | None = None) -> int:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    import MDAnalysis as mda
    import mdtraj as md

    from nnrti.analysis.metrics import pocket_volume_proxy_from_universe

    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--manifest", type=Path, default=Path("manifests/md_manifest.csv"))
    p.add_argument("--window-ns", type=float, default=50.0)
    p.add_argument("--min-ns", type=float, default=80.0)
    p.add_argument("--mutations", nargs="*", default=None)
    p.add_argument("--output", type=Path,
                   default=REPO_ROOT / "results/analysis/pocket_surrogates/tables/"
                                       "nnibp_volume_per_frame.csv")
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
            if traj.unitcell_lengths is None:
                raise ValueError("no box information; minimum image is undefined")
            box = traj.unitcell_lengths * 10.0
            xyz = traj.xyz * 10.0

            u = mda.Universe(str(topo))
            protein = u.select_atoms("protein")
            top = traj.topology
            resids = {r + OFFSET for r in NNIBP_P66}
            pocket_ca = [a.index for a in top.atoms
                         if a.name == "CA" and a.residue.resSeq in resids]
            if not pocket_ca:
                raise ValueError("no pocket CA atoms found")
            prot_idx = [a.index for a in top.atoms if a.residue.name not in ("HOH", "WAT")
                        and a.residue.is_protein]
            if len(prot_idx) != protein.n_atoms:
                raise ValueError(f"protein selections disagree: mdtraj {len(prot_idx)} "
                                 f"vs MDAnalysis {protein.n_atoms}")

            mask = _window_mask(traj.n_frames, total_ns, args.window_ns)
            for f in np.flatnonzero(mask):
                f = int(f)
                frame = xyz[f]
                b = box[f]
                # pocket centroid, imaged against its own first CA
                ref = frame[pocket_ca[0]]
                centre = ref + _mic(frame[pocket_ca] - ref, b).mean(axis=0)
                # every protein atom placed in the image nearest that centroid
                protein.positions = centre + _mic(frame[prot_idx] - centre, b)
                rows.append(dict(
                    mutation=mutation, state=state, replicate=rep, frame=f,
                    run_ns=round(total_ns, 1),
                    pocket_volume=float(pocket_volume_proxy_from_universe(u, "")),
                ))
            logging.info("%s %s rep%02d: %d of %d frames in the last %.0f ns",
                         mutation, state, rep, int(mask.sum()), traj.n_frames,
                         args.window_ns)
        except Exception as exc:
            logging.error("FAILED %s %s rep%s: %s", mutation, state, rep, exc)

    if not rows:
        return 1
    df = pd.DataFrame(rows)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.output, index=False)
    logging.info("wrote %s (%d rows)", args.output, len(df))
    print(df.groupby(["state"]).pocket_volume.describe().to_string())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
