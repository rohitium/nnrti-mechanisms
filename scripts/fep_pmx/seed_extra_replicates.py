#!/usr/bin/env python3
"""Seed additional FEP replicates from the final frame of an existing MD run.

Why
---
FEP replicate error is between-replicate: SEM = sigma_DDG / sqrt(n), so extra
replicates are the one lever that always works. But the FEP leg resolver maps
replicate n -> ``results/md_runs/<label>/rep_{n:02d}/assets/*_start.pdb`` and only
THREE MD replicates exist. There is no rep_04 to point at.

This script creates one by copying a run's ``*_md_final.pdb`` (the last production
frame, ~100 ns from its own start) into the ``rep_{n}/assets/`` path convention.
No format conversion is involved: OpenMM writes ``_md_final.pdb`` with the same
writer, atom count and naming as ``_start.pdb`` -- they differ only in coordinates.

Each leg-replicate needs FOUR structures (source + endpoint) x (holo + apo), so a
leg is only seedable when all four exist. Legs whose apo MD never ran (y181c,
v106i, v106a_l234i as of 2026-08-15 -- their apo run dirs hold only ``assets/``)
cannot be seeded this way; see RUNBOOK / STATUS for the alternatives considered.

Validation (all must pass before anything is written)
-----------------------------------------------------
1. atom count identical to the leg's existing rep_01 start PDB
2. protein is not split across the periodic boundary (extent < 80% of box edge)
3. holo only: ligand present with the expected atom count, and in contact with
   the protein (min heavy-atom distance < 5 A) -- i.e. DOR did not unbind

(3) matters because ``build_solvated_system`` reads the DOR pose directly out of
this PDB (``extract_ligand_coords_nm``); an unbound or drifted ligand would be
baked into the alchemical system silently.

Usage
-----
    python3 scripts/fep_pmx/seed_extra_replicates.py --legs wt_to_V106A            # dry run
    python3 scripts/fep_pmx/seed_extra_replicates.py --legs wt_to_V106A --apply
    python3 scripts/fep_pmx/seed_extra_replicates.py --legs A B --source-rep 2 --dest-rep 5 --apply

Then prepare the new replicate as usual (note REPLICATES is a LIST here):

    LEGS="wt_to_V106A" REPLICATES="4" bash scripts/fep_pmx/prepare_p0_hybrids.sh
    LEGS="wt_to_V106A" REPLICATES="4" bash scripts/fep_pmx/build_p0_systems.sh
    python3 scripts/fep_pmx/prepare_neq.py --legs wt_to_V106A --rep-start 4 --replicates 4 \
        --n-snapshots 100 --force --panel-manifest results/analysis/fep_pmx/neq_<batch>_manifest.csv
"""

from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(REPO_ROOT / "scripts" / "fep_jorgensen"))

LIGAND_RESNAME = "2KW"
SKIP_RESNAMES = {"HOH", "WAT", "SOL", "NA", "CL", "K", "MG", "ZN"}
# PBC-split detection. NOT an absolute fraction of the box: HIV RT is elongated and
# its genuine long axis is ~107 A in a ~128 A box (84%), so any fixed fraction below
# ~0.9 rejects known-good structures -- including the reference _start.pdb itself.
# A real split wraps the protein across the boundary and inflates the extent toward
# the full box edge, so compare against the reference extent, which is known good.
SPLIT_GROWTH_TOLERANCE = 1.25  # reject if any axis exceeds reference extent by >25%
SPLIT_BOX_FRACTION = 0.95      # ...or approaches the box edge outright
CONTACT_CUTOFF_A = 5.0         # ligand must be within this of the protein


def _load(pdb: Path) -> dict:
    """Return atom counts, protein/ligand coords, and box from a PDB."""
    import numpy as np

    protein: list[list[float]] = []
    ligand: list[list[float]] = []
    box: list[float] | None = None
    n_atoms = 0
    for line in pdb.read_text().splitlines():
        if line.startswith("CRYST1"):
            box = [float(line[6:15]), float(line[15:24]), float(line[24:33])]
            continue
        if not line.startswith(("ATOM", "HETATM")):
            continue
        n_atoms += 1
        resname = line[17:20].strip()
        if resname in SKIP_RESNAMES:
            continue
        xyz = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
        (ligand if resname == LIGAND_RESNAME else protein).append(xyz)
    return {
        "n_atoms": n_atoms,
        "protein": np.array(protein),
        "ligand": np.array(ligand),
        "box": box,
        "path": pdb,
    }


def validate(candidate: Path, reference: Path, *, holo: bool) -> list[str]:
    """Return a list of problems; empty means the candidate is usable."""
    import numpy as np

    problems: list[str] = []
    cand = _load(candidate)
    ref = _load(reference)

    if cand["n_atoms"] != ref["n_atoms"]:
        problems.append(f"atom count {cand['n_atoms']} != reference {ref['n_atoms']}")

    if cand["box"] is None:
        problems.append("no CRYST1 box record")
    elif len(cand["protein"]) and len(ref["protein"]):
        extent = cand["protein"].max(0) - cand["protein"].min(0)
        ref_extent = ref["protein"].max(0) - ref["protein"].min(0)
        grown = extent > ref_extent * SPLIT_GROWTH_TOLERANCE
        near_box = extent > np.array(cand["box"]) * SPLIT_BOX_FRACTION
        if (grown | near_box).any():
            problems.append(
                f"protein looks PBC-split: extent {np.round(extent, 1).tolist()} "
                f"vs reference {np.round(ref_extent, 1).tolist()} in box {cand['box']}"
            )

    if holo:
        if not len(cand["ligand"]):
            problems.append(f"no {LIGAND_RESNAME} atoms found")
        else:
            if len(cand["ligand"]) != len(ref["ligand"]):
                problems.append(
                    f"{LIGAND_RESNAME} atom count {len(cand['ligand'])} != reference {len(ref['ligand'])}"
                )
            gap = float(
                np.min(
                    np.linalg.norm(
                        cand["protein"][None, :, :] - cand["ligand"][:, None, :], axis=2
                    )
                )
            )
            if gap > CONTACT_CUTOFF_A:
                problems.append(f"ligand not in contact with protein (min gap {gap:.2f} A)")
    return problems


def _final_pdb_for(start_pdb: Path) -> Path | None:
    """Locate the *_md_final.pdb sitting alongside a rep's assets/ dir."""
    run_dir = start_pdb.parent.parent
    matches = sorted(run_dir.glob("*_md_final.pdb"))
    return matches[0] if matches else None


def _targets_for_leg(leg, source_rep: int, dest_rep: int) -> list[dict]:
    """The four (source/endpoint x holo/apo) copies a leg-replicate needs."""
    return [
        {"what": "holo source", "holo": True,
         "src_start": leg.input_complex_pdb(source_rep), "dst_start": leg.input_complex_pdb(dest_rep)},
        {"what": "holo endpoint", "holo": True,
         "src_start": leg.endpoint_complex_pdb(source_rep), "dst_start": leg.endpoint_complex_pdb(dest_rep)},
        {"what": "apo source", "holo": False,
         "src_start": leg.input_apo_pdb(source_rep), "dst_start": leg.input_apo_pdb(dest_rep)},
        {"what": "apo endpoint", "holo": False,
         "src_start": leg.endpoint_apo_pdb(source_rep), "dst_start": leg.endpoint_apo_pdb(dest_rep)},
    ]


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--legs", nargs="+", required=True, help="leg ids, e.g. wt_to_V106A")
    parser.add_argument("--source-rep", type=int, default=1,
                        help="existing replicate whose final frame is used (default 1)")
    parser.add_argument("--dest-rep", type=int, default=4,
                        help="new replicate index to create (default 4)")
    parser.add_argument("--apply", action="store_true",
                        help="actually write the files (default: dry run)")
    parser.add_argument("--overwrite", action="store_true",
                        help="replace an existing destination start PDB")
    args = parser.parse_args(argv)

    import mutations as m

    plans = m.MANUSCRIPT_PLANS
    items = plans.items() if isinstance(plans, dict) else [(p.genotype, p) for p in plans]
    by_id = {leg.leg_id: leg for _, plan in items for leg in plan.legs}

    unknown = [lid for lid in args.legs if lid not in by_id]
    if unknown:
        parser.error(f"unknown leg id(s): {', '.join(unknown)}")

    print(f"source rep {args.source_rep} -> dest rep {args.dest_rep}   "
          f"({'APPLY' if args.apply else 'DRY RUN'})\n")

    planned: list[tuple[Path, Path]] = []
    blocked = 0
    for lid in args.legs:
        leg = by_id[lid]
        print(f"=== {lid}")
        leg_problems = 0
        leg_actions: list[tuple[Path, Path]] = []
        for t in _targets_for_leg(leg, args.source_rep, args.dest_rep):
            src_start, dst_start = t["src_start"], t["dst_start"]
            final = _final_pdb_for(src_start) if src_start.is_file() else None
            if final is None:
                print(f"  {t['what']:14s} BLOCKED  no *_md_final.pdb beside {src_start.parent.parent}")
                leg_problems += 1
                continue
            if dst_start.is_file() and not args.overwrite:
                print(f"  {t['what']:14s} SKIP     {dst_start} exists (use --overwrite)")
                continue
            problems = validate(final, src_start, holo=t["holo"])
            if problems:
                for p in problems:
                    print(f"  {t['what']:14s} REJECT   {final.name}: {p}")
                leg_problems += 1
                continue
            print(f"  {t['what']:14s} ok       {final.name} -> {dst_start}")
            leg_actions.append((final, dst_start))
        if leg_problems:
            print(f"  -> {lid} NOT seedable ({leg_problems} blocked/rejected); no files written for it\n")
            blocked += 1
            continue
        planned.extend(leg_actions)
        print(f"  -> {lid} seedable ({len(leg_actions)} file(s))\n")

    if not args.apply:
        print(f"DRY RUN: {len(planned)} file(s) would be written; {blocked} leg(s) not seedable.")
        print("Re-run with --apply to write them.")
        return 0

    for final, dst_start in planned:
        dst_start.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(final, dst_start)
        print(f"wrote {dst_start}")
    print(f"\nDone: {len(planned)} file(s) written; {blocked} leg(s) not seedable.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
