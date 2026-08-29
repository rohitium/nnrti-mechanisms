#!/usr/bin/env python3
"""Gate the widened MM/GBSA sampling window on DOR staying in one bound pose.

Why this exists
---------------
MM/GBSA is an ensemble average, so the snapshots should span the whole
post-equilibration trajectory rather than a single terminal slice. That is only
legitimate while the ligand occupies **one** binding mode across the window: if
DOR left its crystallographic pose partway through a run, averaging over the
window would mix two states and the mean would describe neither.

This script measures, over exactly the frames the widened protocol will score,
how far DOR moves. For each run it reports

* ``dor_rmsd_mean/max`` -- DOR heavy-atom RMSD to the crystallographic pose after
  superposing the NNIBP pocket Calphas, in angstrom;
* ``dor_rmsd_self_sd`` -- RMSD of each frame's DOR to the run's own mean DOR
  pose, i.e. the spread of the sampled ensemble independent of any crystal
  reference;
* ``com_sd`` -- spread of the RT-DOR centre-of-mass distance.

A run passes when the ligand stays in one basin: ``dor_rmsd_self_sd`` below
``--max-self-sd`` and ``dor_rmsd_max`` below ``--max-rmsd``. Failures are
reported, never silently dropped -- decide per run whether to narrow its window
or exclude it.

Frame selection mirrors ``compute_mmgbsa_safe --frame-sampling even`` exactly:
the contact-screen whitelist, restricted to frames at or after
``discard_fraction``, sampled evenly to ``--snapshots``.
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import numpy as np
import pandas as pd

#: NNIBP pocket residues used for the superposition, p66 topology numbering
#: (canonical - 3). Source: Cilento, Kirby & Sarafianos, Chem. Rev. 2021.
NNBP_CA_RESIDS_P66 = (97, 98, 100, 103, 104, 105, 176, 178, 185, 186, 187, 224, 226, 231, 315)

LIGAND_RESNAME = "2KW"


def _kabsch_rmsd(mobile: np.ndarray, ref: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return (rotation, translation) superposing ``mobile`` onto ``ref``."""
    mc, rc = mobile.mean(axis=0), ref.mean(axis=0)
    cov = (mobile - mc).T @ (ref - rc)
    v, _s, wt = np.linalg.svd(cov)
    d = np.sign(np.linalg.det(v @ wt))
    rot = v @ np.diag([1.0, 1.0, d]) @ wt
    return rot, rc - mc @ rot


def _even_frames(allowed: np.ndarray, n_frames: int, discard: float, take: int) -> np.ndarray:
    start = int(np.floor(max(0.0, min(0.95, discard)) * n_frames))
    pool = allowed[allowed >= start]
    if pool.size == 0:
        pool = allowed
    take = min(take, pool.size)
    return pool[np.linspace(0, pool.size - 1, num=take, dtype=int)]


def _run_one(row: pd.Series, allowed: np.ndarray, ref_pdb: Path, args) -> dict:
    import MDAnalysis as mda
    from MDAnalysis.analysis import align  # noqa: F401  (import validates install)

    u = mda.Universe(str(row["analysis_topology_pdb"]), str(row["analysis_dcd"]))
    n_frames = len(u.trajectory)
    frames = _even_frames(allowed, n_frames, args.discard_fraction, args.snapshots)

    pocket = u.select_atoms(
        "name CA and segid A and resid " + " ".join(str(r) for r in NNBP_CA_RESIDS_P66)
    )
    if pocket.n_atoms == 0:  # segid naming varies; fall back to the first chain
        pocket = u.select_atoms("name CA and resid " + " ".join(str(r) for r in NNBP_CA_RESIDS_P66))
    lig = u.select_atoms(f"resname {LIGAND_RESNAME} and not name H*")
    if lig.n_atoms == 0 or pocket.n_atoms < 3:
        raise ValueError(f"pocket={pocket.n_atoms} ligand={lig.n_atoms} atoms selected")

    ref = mda.Universe(str(ref_pdb))
    ref_pocket = ref.select_atoms(
        "name CA and segid A and resid " + " ".join(str(r) for r in NNBP_CA_RESIDS_P66)
    )
    ref_lig = ref.select_atoms(f"resname {LIGAND_RESNAME} and not name H*")
    have_ref = ref_pocket.n_atoms == pocket.n_atoms and ref_lig.n_atoms == lig.n_atoms

    poses, rmsd_to_ref, coms = [], [], []
    for f in frames:
        u.trajectory[int(f)]
        rot, trans = _kabsch_rmsd(pocket.positions.copy(), ref_pocket.positions.copy()) \
            if have_ref else (np.eye(3), np.zeros(3))
        lp = lig.positions @ rot + trans
        poses.append(lp)
        if have_ref:
            rmsd_to_ref.append(float(np.sqrt(((lp - ref_lig.positions) ** 2).sum(axis=1).mean())))
        coms.append(float(np.linalg.norm(lig.center_of_mass() - pocket.center_of_mass())))

    poses_arr = np.asarray(poses)
    mean_pose = poses_arr.mean(axis=0)
    self_rmsd = np.sqrt(((poses_arr - mean_pose) ** 2).sum(axis=2).mean(axis=1))

    return {
        "mutation": row["mutation"],
        "replicate": int(row["replicate"]),
        "n_frames_total": n_frames,
        "n_clean": int(allowed.size),
        "n_scored": int(frames.size),
        "first_frame": int(frames[0]),
        "last_frame": int(frames[-1]),
        "dor_rmsd_mean": float(np.mean(rmsd_to_ref)) if rmsd_to_ref else np.nan,
        "dor_rmsd_max": float(np.max(rmsd_to_ref)) if rmsd_to_ref else np.nan,
        "dor_rmsd_self_sd": float(self_rmsd.std(ddof=1)),
        "dor_rmsd_self_max": float(self_rmsd.max()),
        "com_mean": float(np.mean(coms)),
        "com_sd": float(np.std(coms, ddof=1)),
    }


def main() -> int:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    root = Path(__file__).resolve().parents[3]
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--manifest", type=Path, default=root / "manifests/md_manifest.csv")
    p.add_argument("--results-dir", type=Path, default=root / "results")
    p.add_argument("--contact-screen-csv", type=Path,
                   default=root / "results/analysis/binding_energy/frame_contact_screen.csv")
    p.add_argument("--reference-pdb", type=Path,
                   default=root / "results/md_runs/wt/rep_01/wt_rep01_analysis_topology.pdb",
                   help="Pose reference; the minimized WT analysis topology by default.")
    p.add_argument("--snapshots", type=int, default=100)
    p.add_argument("--discard-fraction", type=float, default=0.25)
    p.add_argument("--max-self-sd", type=float, default=1.0,
                   help="Fail a run whose DOR pose spread about its own mean exceeds this (A).")
    p.add_argument("--max-rmsd", type=float, default=3.0,
                   help="Fail a run whose worst DOR RMSD to the reference pose exceeds this (A).")
    p.add_argument("--output", type=Path,
                   default=root / "results/analysis/binding_energy/pose_stability_even_window.csv")
    args = p.parse_args()

    import sys
    sys.path.insert(0, str(root))
    from src.analysis.result_collector import collect_md_results

    md = collect_md_results(args.manifest, args.results_dir)
    screen = pd.read_csv(args.contact_screen_csv)
    screen = screen[(screen["status"] == "ok") & screen["is_clean"].astype(bool)]

    rows = []
    for _, row in md.iterrows():
        mut, rep = str(row["mutation"]), int(row["replicate"])
        allowed = np.array(sorted(
            screen[(screen.mutation == mut) & (screen.replicate == rep)].frame.astype(int)
        ), dtype=int)
        if allowed.size == 0:
            logging.warning("%s rep%d: no clean frames", mut, rep)
            continue
        try:
            rows.append(_run_one(row, allowed, args.reference_pdb, args))
            logging.info("%s rep%d ok", mut, rep)
        except Exception as exc:
            logging.error("%s rep%d failed: %s", mut, rep, exc)

    out = pd.DataFrame(rows).sort_values(["mutation", "replicate"])
    out["pass"] = (out.dor_rmsd_self_sd <= args.max_self_sd) & (
        out.dor_rmsd_max.isna() | (out.dor_rmsd_max <= args.max_rmsd)
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.output, index=False)

    failed = out[~out["pass"]]
    logging.info("wrote %s (%d runs, %d failed)", args.output, len(out), len(failed))
    if not failed.empty:
        logging.warning("runs failing the single-pose gate:\n%s",
                        failed[["mutation", "replicate", "dor_rmsd_self_sd", "dor_rmsd_max"]]
                        .to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
