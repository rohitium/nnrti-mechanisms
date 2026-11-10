#!/usr/bin/env python3
"""Regional Ca RMSF and strand content in 10 ns blocks, apo and holo.

Figure 5 and Supplementary Table 4 carry one RMSF per replicate, averaged over
the whole 100 ns. That is the right unit for the heatmap but it cannot carry a
distribution: three numbers do not support a density, and a violin drawn over
them shows the spread between replicates, not the spread of the quantity.

An RMSF has no per-frame value -- it is a fluctuation about a mean structure,
so it only exists over a window. Blocking the trajectory gives it a genuine
sampling distribution: ten 10 ns blocks per replicate, thirty per genotype and
state.

Block RMSF is systematically SMALLER than the 100 ns RMSF of the same
replicate, because a 10 ns window samples less of the motion. The two are
different quantities and the figure labels them as such; the block
distribution is not an error bar on the published value.

The time axis is state-dependent and this is the one thing in here that must
not be guessed:

  holo  uniform. Every holo run JSON carries analysis_n_frames and
        analysis_dt_ps, repaired on 2026-08-13 against the DCD fingerprint,
        and their product is 100 ns to better than 0.6%. Frame counts run
        139-540, so compare_apo_ensembles.frame_times_ns -- which assumes the
        apo schedule of 200 frames at 0.05 ns followed by 0.5 ns frames --
        reconstructs 160 ns for a 500-frame holo run and is simply the wrong
        axis here. It raises on wt rep_01 rather than returning that silently,
        which is how this was found.

  apo   two strides, via compare_apo_ensembles.frame_times_ns, which is
        checked against the step record.

Neither axis is taken from a DCD header or from ts.time.

The preprocessing is copied from compute_protein_level_changes so the
full-trajectory values this script also writes can be checked against the
published per-residue tables: re-image anchored on p66, drop frames whose
pocket Ca radius of gyration exceeds POCKET_RG_MAX_A, superpose on the pocket
Ca atoms. Superposition for a block is onto that block's own first frame, so a
block RMSF measures fluctuation inside the block rather than drift across the
run.

Usage
-----
    PYTHONPATH=src python -m nnrti.cli.compute_block_metrics
    PYTHONPATH=src python -m nnrti.cli.compute_block_metrics --states holo --genotypes wt
"""
from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path

import numpy as np
import pandas as pd

from nnrti.cli.build_structural_metrics_table import POCKET_RG_MAX_A, REGIONS
from nnrti.cli.compare_apo_ensembles import frame_times_ns
from nnrti.cli.compute_protein_level_changes import (
    MIN_CLEAN_FRAMES,
    NNIBP_P66,
    RESID_OFFSET,
)
from nnrti.cli.find_pocket_surrogates import REPO_ROOT, _run_ns

#: A block with fewer frames than this cannot carry an RMSF worth plotting.
MIN_BLOCK_FRAMES = 8

#: Tolerance on the holo JSON axis, as a fraction of the nominal run length.
HOLO_AXIS_TOL = 0.02


def _display(label: str) -> str:
    return "WT" if label.lower() == "wt" else label.upper().replace("_", "+")


def holo_frame_times_ns(rep_dir: Path, n_frames: int) -> np.ndarray:
    """Per-frame time for a holo replicate, from its own repaired JSON.

    The holo DCDs are uniform in time with a per-replicate stride, so the axis
    is n * dt. The JSON is cross-checked against the step record before it is
    believed, because CLAUDE.md is right that the run JSON is metadata rather
    than a measurement.
    """
    js = next((q for q in sorted(rep_dir.glob("*.json"))), None)
    steps_ns = _run_ns(rep_dir)
    if js is not None:
        meta = json.loads(js.read_text())
        n_json = meta.get("analysis_n_frames")
        dt_ps = meta.get("analysis_dt_ps")
        if n_json and dt_ps:
            if int(n_json) != int(n_frames):
                raise ValueError(
                    f"{rep_dir}: JSON says {n_json} analysis frames but the DCD "
                    f"holds {n_frames}. The trajectory has changed since the "
                    "timing repair; re-run the repair rather than trusting either."
                )
            span = n_frames * float(dt_ps) / 1000.0
            if steps_ns > 0 and span > steps_ns * (1.0 + HOLO_AXIS_TOL):
                raise ValueError(
                    f"{rep_dir}: JSON axis spans {span:.2f} ns, longer than the "
                    f"{steps_ns:.2f} ns the step record records. One of the two is "
                    "wrong; fix the timing repair rather than plotting this."
                )
            if steps_ns > 0 and span < steps_ns * (1.0 - HOLO_AXIS_TOL):
                # G190E rep_03 is the one case: 139 frames at the same 558.659 ps
                # stride as its two sibling replicates, so the analysis DCD stops
                # at 77.6 ns although the run itself reached 50,000,000 steps.
                # The DCD is truncated, not mis-strided. Stretching those frames
                # across 100 ns -- which is what the n_frames <= 200 branch of
                # frame_times_ns does -- would place each one about 29% late and
                # put block 8's frames in block 10. Values averaged over all
                # frames do not depend on the axis, so the published tables are
                # unaffected; only the blocking here is.
                logging.warning(
                    "%s: analysis DCD covers %.1f ns of a %.1f ns run (%d frames "
                    "at %.1f ps); treating it as truncated, not re-strided.",
                    rep_dir, span, steps_ns, n_frames, float(dt_ps))
            return np.arange(1, n_frames + 1) * float(dt_ps) / 1000.0
    if steps_ns <= 0:
        raise ValueError(f"{rep_dir}: no usable time axis")
    return np.arange(1, n_frames + 1) * (steps_ns / n_frames)


def frame_times(state: str, rep_dir: Path, n_frames: int) -> np.ndarray:
    if state == "holo":
        return holo_frame_times_ns(rep_dir, n_frames)
    return frame_times_ns(n_frames, _run_ns(rep_dir))


def _regional(values: np.ndarray, canonical: np.ndarray, chain: np.ndarray,
              prefix: str) -> dict[str, float]:
    """Mean of a per-Ca quantity over each p66 region."""
    out = {}
    for region, resids in REGIONS.items():
        sel = (chain == 0) & np.isin(canonical, list(resids))
        out[f"{prefix}_{region}"] = float(values[sel].mean()) if sel.any() else np.nan
    return out


def _rmsf(xyz: np.ndarray) -> np.ndarray:
    """Per-atom RMSF in Angstrom for coordinates already superposed, (n, a, 3)."""
    return np.sqrt(((xyz - xyz.mean(0)) ** 2).sum(-1).mean(0))


def main(argv: list[str] | None = None) -> int:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    import mdtraj as md

    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--block-ns", type=float, default=10.0)
    p.add_argument("--states", nargs="*", default=["holo", "apo"])
    p.add_argument("--genotypes", nargs="*", default=None)
    p.add_argument("--validate-max-frames", type=int, default=200,
                   help="Frame cap for the full-trajectory validation values, "
                        "matching compute_protein_level_changes.")
    p.add_argument("--output-dir", type=Path,
                   default=REPO_ROOT / "results/analysis/block_metrics/tables")
    args = p.parse_args(argv)

    pocket_resids = [r + RESID_OFFSET for r in NNIBP_P66]
    block_rows, full_rows = [], []

    for state in args.states:
        root = Path("results/md_runs/apo") if state == "apo" else Path("results/md_runs")
        dirs = sorted(d for d in root.glob("*/rep_*")
                      if (state == "apo" or "apo" not in d.parts)
                      and list(d.glob("*analysis.dcd")))
        if args.genotypes:
            want = {g.lower() for g in args.genotypes}
            dirs = [d for d in dirs if d.parent.name.lower() in want]
        for d in dirs:
            label, rep = d.parent.name, int(d.name.split("_")[-1])
            try:
                dcd = sorted(d.glob("*analysis.dcd"))[0]
                top = sorted(d.glob("*analysis_topology.pdb"))[0]
                t = md.load(str(dcd), top=str(top))
                t.image_molecules(inplace=True,
                                  anchor_molecules=[set(t.topology.chain(0).atoms)])
                times = frame_times(state, d, t.n_frames)

                pocket_ca = [a.index for a in t.topology.atoms
                             if a.name == "CA" and a.residue.chain.index == 0
                             and a.residue.resSeq in pocket_resids]
                if len(pocket_ca) != len(NNIBP_P66):
                    raise ValueError(f"found {len(pocket_ca)} pocket CA atoms, "
                                     f"expected {len(NNIBP_P66)}")
                pk = t.xyz[:, pocket_ca, :] * 10.0
                rg = np.sqrt(((pk - pk.mean(axis=1, keepdims=True)) ** 2)
                             .sum(-1).mean(-1))
                keep = rg <= POCKET_RG_MAX_A
                if keep.sum() < MIN_CLEAN_FRAMES:
                    raise ValueError(f"only {int(keep.sum())} of {t.n_frames} frames "
                                     "have an intact pocket")
                if not keep.all():
                    logging.info("  %s %s rep%02d: dropped %d of %d frames, pocket "
                                 "Rg > %.0f A", state, label, rep,
                                 int((~keep).sum()), t.n_frames, POCKET_RG_MAX_A)
                t = t[keep]
                times = times[keep]

                ca = [a for a in t.topology.atoms if a.name == "CA"]
                ca_idx = [a.index for a in ca]
                canonical = np.array([a.residue.resSeq - RESID_OFFSET for a in ca])
                chain = np.array([a.residue.chain.index for a in ca])
                # compute_dssp returns one column per residue in the topology,
                # the ligand included as 'NA', so the protein columns are taken
                # by residue index rather than by position in a filtered list.
                dssp_all = md.compute_dssp(t, simplified=True)
                prot = [r for r in t.topology.residues if r.is_protein]
                dssp = dssp_all[:, [r.index for r in prot]]
                ss_canonical = np.array([r.resSeq - RESID_OFFSET for r in prot])
                ss_chain = np.array([r.chain.index for r in prot])

                # ---- full trajectory, the published protocol ----------------
                sub = (np.linspace(0, t.n_frames - 1, args.validate_max_frames).astype(int)
                       if t.n_frames > args.validate_max_frames
                       else np.arange(t.n_frames))
                tf = t[sub]
                tf.superpose(tf, 0, atom_indices=pocket_ca)
                row = {"state": state, "genotype": _display(label), "replicate": rep,
                       "n_frames": len(sub), "total_ns": round(float(times[-1]), 2)}
                row.update(_regional(_rmsf(tf.xyz[:, ca_idx, :] * 10.0),
                                     canonical, chain, "rmsf_A"))
                row.update(_regional((dssp[sub] == "E").mean(0),
                                     ss_canonical, ss_chain, "strand_frac"))
                full_rows.append(row)

                # ---- 10 ns blocks ------------------------------------------
                edges = np.arange(0.0, times[-1] + args.block_ns, args.block_ns)
                for b in range(len(edges) - 1):
                    lo, hi = edges[b], edges[b + 1]
                    sel = np.flatnonzero((times > lo) & (times <= hi))
                    if len(sel) < MIN_BLOCK_FRAMES:
                        continue
                    tb = t[sel]
                    tb.superpose(tb, 0, atom_indices=pocket_ca)
                    brow = {"state": state, "genotype": _display(label),
                            "replicate": rep, "block": b + 1,
                            "t_start_ns": round(float(lo), 2),
                            "t_end_ns": round(float(hi), 2), "n_frames": len(sel)}
                    brow.update(_regional(_rmsf(tb.xyz[:, ca_idx, :] * 10.0),
                                          canonical, chain, "rmsf_A"))
                    brow.update(_regional((dssp[sel] == "E").mean(0),
                                          ss_canonical, ss_chain, "strand_frac"))
                    block_rows.append(brow)
                logging.info("%s %-12s rep%02d  %d frames over %.1f ns  %d blocks",
                             state, label, rep, t.n_frames, times[-1],
                             sum(1 for r in block_rows
                                 if r["genotype"] == _display(label)
                                 and r["replicate"] == rep and r["state"] == state))
            except Exception as exc:
                logging.error("FAILED %s %s rep%02d: %s", state, label, rep, exc)

    if not block_rows:
        return 1
    args.output_dir.mkdir(parents=True, exist_ok=True)
    blocks = pd.DataFrame(block_rows).round(4)
    blocks.to_csv(args.output_dir / "block_metrics.csv", index=False)
    full = pd.DataFrame(full_rows).round(4)
    full.to_csv(args.output_dir / "full_trajectory_metrics.csv", index=False)
    logging.info("wrote %s (%d blocks, %d replicates)", args.output_dir,
                 len(blocks), len(full))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
