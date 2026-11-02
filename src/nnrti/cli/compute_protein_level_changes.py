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

#: Frames whose pocket Ca radius of gyration exceeds this are dropped. The
#: RMSF below is measured after superposing on those same pocket Ca atoms, so a
#: pocket that has come apart destroys the reference frame and every residue in
#: the enzyme then appears to fly. apo g190e rep_02 is the case in point: its
#: pocket Rg runs to 42 A and all 955 residues read above 6 A RMSF. The
#: threshold is the one compare_apo_ensembles.py already uses, so the two
#: analyses exclude the same frames.
POCKET_RG_MAX_A = 10.0

#: A replicate with fewer surviving frames than this is reported rather than
#: scored, so a mostly-unfolded run cannot contribute a number that looks fine.
MIN_CLEAN_FRAMES = 20

NNIBP_P66 = (100, 101, 103, 106, 107, 108, 179, 181, 188, 189, 190, 227, 229, 234, 318)

#: p66 subdomain boundaries, canonical numbering, as used in
#: compare_dor_crystal_poses.py. Kept identical so the two analyses partition
#: the enzyme the same way.
_FINGERS = tuple(range(1, 86)) + tuple(range(118, 156))
_PALM = tuple(range(86, 118)) + tuple(range(156, 237))
_THUMB = tuple(range(237, 319))
_CONNECTION = tuple(range(319, 427))

REGIONS = {
    "NNIBP lining": NNIBP_P66,
    "primer grip 225-235": tuple(range(225, 236)),
    "beta6 100-110": tuple(range(100, 111)),
    "beta9-beta10 179-190": tuple(range(179, 191)),
    # The catalytic motif and the subdomains whose relative motion defines the
    # open/closed transition: a mutation can leave the pocket lining alone and
    # still change how the enzyme breathes.
    "YMDD 183-186": (183, 184, 185, 186),
    "p66 fingers": _FINGERS,
    "p66 palm": _PALM,
    "p66 thumb": _THUMB,
    "p66 connection": _CONNECTION,
}


def main(argv: list[str] | None = None) -> int:
    import mdtraj as md

    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--state", choices=("holo", "apo"), default="holo",
                   help="Which leg to analyse. Sets the default runs-root and "
                        "output-dir, and is recorded in a structure_type column "
                        "so the two legs can be concatenated unambiguously.")
    p.add_argument("--runs-root", type=Path, default=None)
    p.add_argument("--output-dir", type=Path, default=None)
    p.add_argument("--max-frames", type=int, default=200,
                   help="Subsample to at most this many frames per replicate.")
    args = p.parse_args(argv)

    # Defaults depend on the leg, and the apo runs live one level deeper. The
    # holo filter must still exclude results/md_runs/apo, but that same test
    # would reject every apo directory, so it is applied only for holo.
    if args.runs_root is None:
        args.runs_root = (Path("results/md_runs/apo") if args.state == "apo"
                          else Path("results/md_runs"))
    if args.output_dir is None:
        args.output_dir = Path(f"results/analysis/protein_level"
                               f"{'_apo' if args.state == 'apo' else ''}")

    dirs = sorted(d for d in args.runs_root.glob("*/rep_*")
                  if (args.state == "apo" or "apo" not in d.parts)
                  and (list(d.glob("*analysis.dcd"))))
    if not dirs:
        raise SystemExit(f"no {args.state} replicates found under {args.runs_root}")
    print(f"{args.state}: {len(dirs)} replicates under {args.runs_root}")

    rmsf_rows, ss_rows = [], []
    for d in dirs:
        label, rep = d.parent.name, int(d.name.split("_")[-1])
        dcd = sorted(d.glob("*analysis.dcd"))[0]
        top = sorted(d.glob("*analysis_topology.pdb"))[0]
        try:
            t = md.load(str(dcd), top=str(top))
            # Re-image in memory, for BOTH legs, anchored on p66.
            #
            # make_molecules_whole is NOT enough here and is in fact a no-op on
            # these trajectories: p66 and p51 are separate molecules, each
            # already whole, that drift into different periodic images. The
            # superposition this RMSF depends on then sees p51 displaced by a
            # box vector. Measured on apo wt rep_02: raw and
            # make_molecules_whole both give 24.1 A mean Ca RMSF with a 41 A
            # maximum, while anchoring on p66 gives 1.59 A / 4.61 A -- matching
            # that run's *_analysis_pbcfix.dcd sidecar exactly.
            #
            # Correcting here rather than reading sidecars matters because the
            # apo runs have sidecars and the holo runs have none, so reading
            # them would treat the two legs differently -- which is precisely
            # the comparison this table exists to make. This writes nothing
            # into the run directories, holo least of all.
            t.image_molecules(inplace=True,
                              anchor_molecules=[set(t.topology.chain(0).atoms)])
        except Exception as e:
            print(f"  {label} rep{rep}: load failed ({e})", flush=True)
            continue
        if t.n_frames > args.max_frames:
            t = t[np.linspace(0, t.n_frames - 1, args.max_frames).astype(int)]

        # Chain 0 ONLY. NNIBP_P66 are p66 residue numbers, but p51 carries the
        # same numbers, so a chain-blind match returns 30 atoms instead of 15 --
        # the 15 pocket residues plus 15 unrelated p51 residues ~30 A away. The
        # superposition anchor was then a p66/p51 composite rather than the
        # pocket, and since the two subunits move relative to one another, that
        # put their relative motion into every residue's RMSF. The giveaway is
        # the radius of gyration: 31.2 A for the chain-blind selection against
        # 7.7 A for the pocket, the latter matching apo_frame_metrics.
        pocket_ca = [a.index for a in t.topology.atoms
                     if a.name == "CA" and a.residue.chain.index == 0
                     and a.residue.resSeq in [r + RESID_OFFSET for r in NNIBP_P66]]
        if len(pocket_ca) != len(NNIBP_P66):
            print(f"  {label} rep{rep}: expected {len(NNIBP_P66)} pocket CA atoms, "
                  f"found {len(pocket_ca)} -- selection is wrong, skipping", flush=True)
            continue
        if not pocket_ca:
            print(f"  {label} rep{rep}: no pocket CA found", flush=True)
            continue
        # Drop frames whose pocket has come apart BEFORE superposing on it.
        pk = t.xyz[:, pocket_ca, :] * 10.0
        rg = np.sqrt(((pk - pk.mean(axis=1, keepdims=True)) ** 2).sum(-1).mean(-1))
        keep = rg <= POCKET_RG_MAX_A
        if keep.sum() < MIN_CLEAN_FRAMES:
            print(f"  {label} rep{rep}: SKIPPED -- only {int(keep.sum())} of "
                  f"{t.n_frames} frames have an intact pocket "
                  f"(Rg {rg.min():.1f}-{rg.max():.1f} A)", flush=True)
            continue
        if not keep.all():
            print(f"  {label} rep{rep}: dropped {int((~keep).sum())} of {t.n_frames} "
                  f"frames with pocket Rg > {POCKET_RG_MAX_A:.0f} A", flush=True)
            t = t[keep]
        t.superpose(t, 0, atom_indices=pocket_ca)

        ca = [a for a in t.topology.atoms if a.name == "CA"]
        idx = [a.index for a in ca]
        xyz = t.xyz[:, idx, :] * 10.0
        rmsf = np.sqrt(((xyz - xyz.mean(0)) ** 2).sum(-1).mean(0))
        for a, v in zip(ca, rmsf):
            rmsf_rows.append({"structure_type": args.state,
                              "safe_label": label, "replicate": rep,
                              "resSeq": a.residue.resSeq,
                              "canonical": a.residue.resSeq - RESID_OFFSET,
                              "chain": a.residue.chain.index, "rmsf": float(v)})

        try:
            dssp = md.compute_dssp(t, simplified=True)
            residues = [r for r in t.topology.residues if r.is_protein]
            for j, r in enumerate(residues):
                col = dssp[:, j]
                ss_rows.append({"structure_type": args.state,
                                "safe_label": label, "replicate": rep,
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
        row = {"structure_type": args.state, "safe_label": lbl}
        for name, resids in REGIONS.items():
            sel = delta[[i for i in delta.index if i[1] in resids and i[0] == 0]]
            row[name] = round(float(sel.mean()), 3) if len(sel) else float("nan")
        rows.append(row)
    summ = pd.DataFrame(rows).set_index(["structure_type", "safe_label"]).sort_index()
    summ.to_csv(out / "rmsf_delta_by_region.csv")

    print("\n=== mean Ca RMSF change vs WT, by region (A) ===")
    print(summ.to_string())
    print(f"\nwrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
