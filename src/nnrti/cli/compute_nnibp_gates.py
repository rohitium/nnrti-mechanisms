#!/usr/bin/env python3
"""Distances across the NNIBP entrance, per frame, in both states.

The literature mechanism for K103N is not a change in cavity volume. It is a
hydrogen bond: the Asn103 side chain reaches Tyr188 and latches the mouth of
the pocket shut, stabilising the closed unliganded form and slowing inhibitor
entry (Hsiou et al.; Ren et al.). A pocket-volume proxy measures the wrong
coordinate for that, which is why K103N looks identical to wild type in
Supplementary Figure 6 panel A while the literature says it should not.

This measures the mouth directly, as minimum heavy-atom distances between the
residue pairs that line it:

    103-188   the K103N latch itself
    101-188   the vertical gate
    106-181   the lateral gate
    100-227   the floor of the entrance channel
    103-225   the rim toward the primer grip

Reported in both states. The unliganded state is the one that bears on entry;
the bound state is reported alongside so a gate that is simply narrow in this
construct can be told from one that closes only when the drug is absent.

Minimum heavy-atom distance rather than a named atom pair, because the
mutations change which atoms exist: K103N replaces the lysine amine with an
amide, Y188L removes the hydroxyl entirely. A fixed atom pair would be
undefined in exactly the genotypes the measurement is for. Pairs involving a
residue that a genotype deletes are reported as missing rather than guessed.

Usage
-----
    PYTHONPATH=src python -m nnrti.cli.compute_nnibp_gates
"""
from __future__ import annotations

import argparse
import logging
from pathlib import Path

import numpy as np
import pandas as pd

from nnrti.cli.compute_block_metrics import _display, frame_times
from nnrti.cli.compute_protein_level_changes import RESID_OFFSET

#: A gate wider than this means the pocket is no longer a pocket. These are
#: side chains that line one mouth of a cavity; the largest value any intact
#: replicate reaches is under 15 A, and the only frames above the cutoff are
#: G190E unliganded replicate 2 from 20 ns on, which is the unfolding already
#: documented for that run (Supplementary Figure 6 caption). Keeping them
#: turns that genotype's channel floor into 30 A and its primer-grip rim into
#: 29 A, both of which are measurements of a dissociated chain.
GATE_MAX_A = 20.0

#: (canonical residue a, canonical residue b, label). p66 both sides.
GATES = [
    (103, 188, "103-188 latch"),
    (101, 188, "101-188 vertical gate"),
    (106, 181, "106-181 lateral gate"),
    (100, 227, "100-227 channel floor"),
    (103, 225, "103-225 primer-grip rim"),
]


def _heavy(top, canonical: int) -> np.ndarray | None:
    """Heavy-atom indices of one p66 residue, or None if it is absent."""
    res = [r for r in top.residues
           if r.resSeq == canonical + RESID_OFFSET and r.chain.index == 0]
    if len(res) != 1:
        return None
    idx = [a.index for a in res[0].atoms if a.element.symbol != "H"]
    return np.asarray(idx, dtype=int) if idx else None


def main(argv: list[str] | None = None) -> int:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    import mdtraj as md

    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--states", nargs="*", default=["holo", "apo"])
    p.add_argument("--genotypes", nargs="*", default=None)
    p.add_argument("--output-dir", type=Path,
                   default=Path("results/analysis/nnibp_gates/tables"))
    args = p.parse_args(argv)

    rows = []
    for state in args.states:
        root = Path("results/md_runs/apo") if state == "apo" else Path("results/md_runs")
        dirs = sorted(d for d in root.glob("*/rep_*")
                      if (state == "apo" or "apo" not in d.parts)
                      and "pose" not in str(d) and "solvent" not in str(d)
                      and "collapse" not in str(d)
                      and list(d.glob("*analysis.dcd")))
        if args.genotypes:
            want = {g.lower() for g in args.genotypes}
            dirs = [d for d in dirs if d.parent.name.lower() in want]
        for d in dirs:
            label, rep = d.parent.name, int(d.name.split("_")[-1])
            try:
                dcd = sorted(q for q in d.glob("*analysis.dcd")
                             if "pbcfix" not in q.name)[0]
                top = sorted(q for q in d.glob("*analysis_topology.pdb")
                             if "from_md" not in q.name)[0]
                t = md.load(str(dcd), top=str(top))
                # A gate is a pair of side chains a few Angstrom apart inside
                # one subunit; make the molecule whole so neither is reported
                # across a box face.
                t.make_molecules_whole(inplace=True)
                times = frame_times(state, d, t.n_frames)
                for a, b, name in GATES:
                    ia, ib = _heavy(t.topology, a), _heavy(t.topology, b)
                    if ia is None or ib is None:
                        continue
                    pairs = np.array([(i, j) for i in ia for j in ib], dtype=int)
                    dist = md.compute_distances(t, pairs, periodic=False)
                    mind = dist.min(axis=1) * 10.0  # nm -> A
                    for f, (v, tt) in enumerate(zip(mind, times)):
                        rows.append(dict(state=state, genotype=_display(label),
                                         replicate=rep, gate=name, frame=f,
                                         time_ns=round(float(tt), 3),
                                         min_dist_A=round(float(v), 3)))
                logging.info("%s %-12s rep%02d  %d frames", state, label, rep,
                             t.n_frames)
            except Exception as exc:
                logging.error("FAILED %s %s rep%02d: %s", state, label, rep, exc)

    if not rows:
        return 1
    df = pd.DataFrame(rows)
    # Drop whole frames, not single gates: if the pocket has come apart, every
    # distance measured on it that frame is meaningless, not just the widest.
    bad = df.loc[df.min_dist_A > GATE_MAX_A, ["state", "genotype", "replicate", "frame"]]
    if len(bad):
        key = set(map(tuple, bad.drop_duplicates().to_numpy()))
        before = len(df)
        mask = [(s_, g, r, f) not in key for s_, g, r, f in
                zip(df.state, df.genotype, df.replicate, df.frame)]
        df = df[mask]
        logging.warning("dropped %d rows in %d frame(s) with a gate wider than "
                        "%.0f A: %s", before - len(df), len(key), GATE_MAX_A,
                        sorted({(g, r) for _s, g, r, _f in key}))
    args.output_dir.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.output_dir / "nnibp_gate_distances_per_frame.csv", index=False)
    per_rep = (df.groupby(["state", "genotype", "replicate", "gate"], as_index=False)
               .min_dist_A.mean())
    per_rep.to_csv(args.output_dir / "nnibp_gate_distances_per_replicate.csv",
                   index=False)
    logging.info("wrote %s (%d frames)", args.output_dir, len(df))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
