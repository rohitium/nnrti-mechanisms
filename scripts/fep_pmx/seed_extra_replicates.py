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


def validate(candidate: Path, reference: Path, *, holo: bool, solute_only: bool = False) -> list[str]:
    """Return a list of problems; empty means the candidate is usable.

    ``solute_only`` marks a frame extracted from the stripped analysis trajectory:
    it legitimately lacks solvent, so total atom count cannot be compared against
    the solvated reference -- the protein atom count is compared instead.
    """
    import numpy as np

    problems: list[str] = []
    cand = _load(candidate)
    ref = _load(reference)

    if solute_only:
        if len(cand["protein"]) != len(ref["protein"]):
            problems.append(
                f"protein atom count {len(cand['protein'])} != reference {len(ref['protein'])}"
            )
    elif cand["n_atoms"] != ref["n_atoms"]:
        problems.append(f"atom count {cand['n_atoms']} != reference {ref['n_atoms']}")

    if cand["box"] is None and not solute_only:
        problems.append("no CRYST1 box record")
    elif len(cand["protein"]) and len(ref["protein"]):
        extent = cand["protein"].max(0) - cand["protein"].min(0)
        ref_extent = ref["protein"].max(0) - ref["protein"].min(0)
        grown = extent > ref_extent * SPLIT_GROWTH_TOLERANCE
        near_box = (
            extent > np.array(cand["box"]) * SPLIT_BOX_FRACTION
            if cand["box"] is not None
            else np.zeros_like(extent, dtype=bool)
        )
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


def _extract_last_frame(start_pdb: Path, out_pdb: Path) -> Path | None:
    """Fall back to the last frame of the analysis trajectory.

    ``_md_final.pdb`` is written only AFTER ``prod.step(remaining_steps)`` returns
    (md_protocol.py), so a job killed by the SLURM wall leaves a checkpoint and an
    incrementally-written ``*_analysis.dcd`` but no final PDB. Since 100 ns takes
    several 12 h segments, that would otherwise block seeding for days.

    The analysis DCD is written with ``atom_indices=solute_idx`` -- solute only
    (protein, plus ligand when holo), which is exactly what the seed needs: the
    downstream ``extract_protein_only`` discards solvent anyway and
    ``build_p0_systems`` re-solvates from scratch.

    PBC is corrected with MDTraj ``make_molecules_whole`` -- the repo's
    authoritative pattern. MDTraj traverses the bond graph, whereas MDAnalysis
    ``unwrap`` fails when the protein COM sits near a box boundary.
    """
    run_dir = start_pdb.parent.parent
    # Prefer an already-PBC-corrected trajectory when one exists.
    traj_path = None
    for pattern in ("*_analysis_pbcfix.dcd", "*_analysis.dcd"):
        matches = sorted(run_dir.glob(pattern))
        if matches:
            traj_path = matches[0]
            break
    tops = sorted(run_dir.glob("*_analysis_topology.pdb"))
    if traj_path is None or not tops:
        return None
    top_path = tops[0]

    import mdtraj as md

    traj = md.load(str(traj_path), top=str(top_path))
    if traj.n_frames == 0:
        return None
    frame = traj[-1]
    frame.make_molecules_whole(inplace=True)
    out_pdb.parent.mkdir(parents=True, exist_ok=True)
    frame.save_pdb(str(out_pdb))
    return out_pdb


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
    parser.add_argument("--no-trajectory", action="store_true",
                        help="require *_md_final.pdb; do not fall back to the last "
                             "analysis-trajectory frame")
    args = parser.parse_args(argv)

    import tempfile
    tmpdir_ctx = tempfile.TemporaryDirectory(prefix="seed_reps_")
    tmpdir = Path(tmpdir_ctx.name)

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
            if dst_start.is_file() and not args.overwrite:
                print(f"  {t['what']:14s} SKIP     {dst_start} exists (use --overwrite)")
                continue

            source = _final_pdb_for(src_start) if src_start.is_file() else None
            solute_only = False
            if source is None and src_start.is_file() and not args.no_trajectory:
                # Run walled out before writing _md_final.pdb: use the last
                # analysis-trajectory frame instead.
                tmp = tmpdir / f"{dst_start.stem}.extracted.pdb"
                source = _extract_last_frame(src_start, tmp)
                solute_only = source is not None
            if source is None:
                print(f"  {t['what']:14s} BLOCKED  no *_md_final.pdb and no usable "
                      f"analysis trajectory under {src_start.parent.parent}")
                leg_problems += 1
                continue

            problems = validate(source, src_start, holo=t["holo"], solute_only=solute_only)
            if problems:
                for p in problems:
                    print(f"  {t['what']:14s} REJECT   {source.name}: {p}")
                leg_problems += 1
                continue
            origin = "traj frame" if solute_only else "md_final"
            print(f"  {t['what']:14s} ok [{origin:9s}] {source.name} -> {dst_start}")
            leg_actions.append((source, dst_start))
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
