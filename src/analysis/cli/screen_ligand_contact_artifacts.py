#!/usr/bin/env python3
"""Flag trajectory frames carrying unphysical ligand-protein heavy-atom contacts.

Some analysis DCDs contain frames in which doravirine and the NNIBP are placed
far too close together -- heavy-atom separations of 0.9-2.4 A, well inside any
chemically possible contact. Both molecules are internally intact in these
frames (bond lengths are normal); only their relative placement is wrong, so the
defect is in the recorded trajectory rather than in the underlying simulation.

Scoring such a frame with MM/GBSA produces an enormous positive r^-12 Lennard-Jones
term that swamps the real interaction energy, which is the dominant source of
replicate-to-replicate scatter in the binding-energy panel.

The per-frame minimum ligand-protein heavy-atom distance is strongly bimodal --
artifact frames sit below ~2.4 A and physical frames above ~2.6 A, with the gap
between essentially unpopulated -- so the default 2.5 A threshold falls in empty
space rather than cutting through a distribution.

Writes one row per frame so downstream tools can choose their own sampling.
"""

from __future__ import annotations

import argparse
import concurrent.futures as cf
import logging
from pathlib import Path

import numpy as np
import pandas as pd

DEFAULT_OUTPUT = Path("results/analysis/binding_energy/frame_contact_screen.csv")


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--manifest", type=Path, default=Path("manifests/md_manifest.csv"))
    parser.add_argument("--output-csv", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--ligand-resname", type=str, default="2KW")
    parser.add_argument(
        "--threshold-a",
        type=float,
        default=2.5,
        help="Frames whose minimum ligand-protein heavy-atom distance falls below this are flagged.",
    )
    parser.add_argument(
        "--neighbor-cutoff-nm",
        type=float,
        default=0.35,
        help="Neighbour-list radius; must exceed --threshold-a or short contacts will be missed.",
    )
    parser.add_argument("--workers", type=int, default=4)
    return parser.parse_args()


def _screen_one(task: dict) -> list[dict]:
    import mdtraj as md

    top, dcd = Path(task["topology"]), Path(task["trajectory"])
    if not (top.exists() and dcd.exists()):
        return [dict(mutation=task["mutation"], replicate=task["replicate"], frame=-1,
                     min_dist_a=np.nan, is_clean=False, status="missing")]
    try:
        traj = md.load(str(dcd), top=str(top))
        traj.make_molecules_whole(inplace=True)
        # Residue names beginning with a digit must be quoted or the selection
        # silently matches nothing.
        ligand = traj.topology.select(f"resname '{task['ligand_resname']}' and element != H")
        protein = traj.topology.select("protein and element != H")
        if ligand.size == 0 or protein.size == 0:
            raise ValueError(f"empty selection: ligand={ligand.size} protein={protein.size}")

        cutoff = float(task["neighbor_cutoff_nm"])
        threshold = float(task["threshold_a"])
        rows = []
        for frame in range(traj.n_frames):
            single = traj.slice([frame])
            near = md.compute_neighbors(single, cutoff, ligand, haystack_indices=protein)[0]
            if near.size == 0:
                # Nothing within the neighbour radius, so nothing near the threshold.
                min_dist = np.nan
            else:
                pairs = np.array(np.meshgrid(ligand, near)).T.reshape(-1, 2)
                min_dist = float(md.compute_distances(single, pairs)[0].min() * 10.0)
            rows.append(dict(mutation=task["mutation"], replicate=task["replicate"], frame=frame,
                             min_dist_a=min_dist,
                             is_clean=bool(np.isnan(min_dist) or min_dist >= threshold),
                             status="ok"))
        return rows
    except Exception as exc:  # noqa: BLE001 - recorded per run, screen continues
        return [dict(mutation=task["mutation"], replicate=task["replicate"], frame=-1,
                     min_dist_a=np.nan, is_clean=False, status=f"error: {exc}")]


def main() -> int:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    args = _parse_args()
    if args.neighbor_cutoff_nm * 10.0 <= args.threshold_a:
        raise ValueError(
            f"--neighbor-cutoff-nm ({args.neighbor_cutoff_nm} nm) must exceed "
            f"--threshold-a ({args.threshold_a} A) or contacts below the threshold are missed"
        )

    manifest = pd.read_csv(args.manifest)
    tasks = []
    for _, row in manifest.iterrows():
        safe, rep = str(row["safe_label"]), int(row["replicate"])
        base = Path(f"results/md_runs/{safe}/rep_{rep:02d}")
        tasks.append(dict(
            mutation=str(row["mutation"]), replicate=rep,
            topology=base / f"{safe.lower()}_rep{rep:02d}_analysis_topology.pdb",
            trajectory=base / f"{safe.lower()}_rep{rep:02d}_analysis.dcd",
            ligand_resname=args.ligand_resname,
            threshold_a=args.threshold_a,
            neighbor_cutoff_nm=args.neighbor_cutoff_nm,
        ))

    rows: list[dict] = []
    with cf.ProcessPoolExecutor(max_workers=max(1, args.workers)) as pool:
        for done, result in enumerate(pool.map(_screen_one, tasks), 1):
            rows.extend(result)
            logging.info("[%d/%d] %s rep%d", done, len(tasks),
                         result[0]["mutation"], result[0]["replicate"])

    frames = pd.DataFrame(rows).sort_values(["mutation", "replicate", "frame"]).reset_index(drop=True)
    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    frames.to_csv(args.output_csv, index=False)

    ok = frames[frames.status == "ok"]
    per_run = ok.groupby(["mutation", "replicate"]).is_clean.agg(["sum", "count"])
    logging.info("Wrote %s", args.output_csv)
    logging.info("runs screened: %d", len(per_run))
    logging.info("frames: %d clean / %d total (%.1f%% flagged)",
                 int(ok.is_clean.sum()), len(ok), 100 * (1 - ok.is_clean.mean()))
    starved = per_run[per_run["sum"] < 20]
    if not starved.empty:
        logging.warning("runs with fewer than 20 clean frames:\n%s", starved.to_string())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
