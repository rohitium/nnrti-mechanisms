"""NNIBP pocket volume per production replicate.

Writes the per-replicate volumes that Table 3 reports, and their per-genotype
means. The volume itself is computed by
:func:`nnrti.analysis.metrics.pocket_volume_proxy_from_universe`: a 0.75 A grid
over a 10 A sphere on the Calpha centroid of the sixteen pocket-lining residues,
counting points further from every protein heavy atom than that atom's van der
Waals radius plus a 1.4 A probe.

Usage:
    python -m nnrti.cli.compute_pocket_volume
    python -m nnrti.cli.compute_pocket_volume --mutations V106A Y188L
"""

from __future__ import annotations

import argparse
import json
import logging
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd

from ..analysis.md_timing import infer_production_ns
from ..analysis.metrics import pocket_volume_proxy_from_universe
from ..analysis.pbc import load_mdtraj_trajectory, pbcfix_dcd_for, raw_analysis_dcd_for
from ..analysis.result_collector import _prepare_profile_jobs, collect_md_results

LOGGER = logging.getLogger("pocket_volume")


def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--manifest", type=Path, default=Path("manifests/md_manifest.csv"))
    p.add_argument("--output-dir", type=Path, default=Path("results/analysis/pocket_volume"))
    p.add_argument("--mutations", nargs="*", default=None)
    p.add_argument("--pocket-max-frames", type=int, default=40, help="Frames scored per replicate")
    p.add_argument("--workers", type=int, default=4)
    return p.parse_args()


def _display_mutation(raw: str) -> str:
    return str(raw).replace("_", "+").strip()


def _frame_subset(n_frames: int, max_frames: int, stride: int) -> np.ndarray:
    base = np.arange(0, n_frames, max(1, stride), dtype=int)
    if max_frames is not None and len(base) > max_frames:
        pick = np.linspace(0, len(base) - 1, num=max_frames, dtype=int)
        base = base[pick]
    return base


def _process_job(job: dict, pocket_max_frames: int) -> dict:
    import MDAnalysis as mda

    traj_path = Path(job["trajectory"])
    topo_path = Path(job["topology"])
    traj = load_mdtraj_trajectory(traj_path, topo_path)
    n_frames = int(traj.n_frames)
    if n_frames < 2:
        raise ValueError("need >=2 frames")

    idx = _frame_subset(n_frames, max_frames=pocket_max_frames,
                        stride=max(1, n_frames // pocket_max_frames))
    u = mda.Universe(str(topo_path), str(traj_path))
    vals = []
    for fi in idx:
        u.trajectory[int(fi)]
        vals.append(float(pocket_volume_proxy_from_universe(u)))
    vals = np.asarray(vals, dtype=float)

    return {
        "mutation": _display_mutation(job["mutation"]),
        "replicate": int(job["replicate"]),
        "pocket_volume_mean": float(np.mean(vals)),
        "pocket_volume_std": float(np.std(vals)),
        "total_ns": float(job["production_ps"] / 1000.0),
    }


def main() -> int:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    args = _parse_args()
    tables = args.output_dir / "tables"
    cfg = args.output_dir / "config"
    for d in (tables, cfg):
        d.mkdir(parents=True, exist_ok=True)

    run_df = collect_md_results(args.manifest)
    if args.mutations:
        keep = {_display_mutation(m) for m in args.mutations}
        run_df = run_df[run_df["mutation"].map(_display_mutation).isin(keep)].copy()

    jobs = []
    for job in _prepare_profile_jobs(run_df):
        traj = Path(job["trajectory"])
        if "pbcfix" not in traj.name:
            pbc = pbcfix_dcd_for(raw_analysis_dcd_for(traj) if traj.exists() else traj)
            if not pbc.exists():
                LOGGER.warning("skip %s rep%s: missing pbcfix", job["mutation"], job["replicate"])
                continue
            job["trajectory"] = str(pbc)
        jobs.append(job)
    if not jobs:
        LOGGER.error("no jobs with pbcfix DCDs")
        return 1

    for job in jobs:
        traj_path = Path(job["trajectory"])
        raw_dcd = raw_analysis_dcd_for(traj_path)
        safe, rep = str(job["safe_label"]), int(job["replicate"])
        call = infer_production_ns(
            dcd_path=raw_dcd if raw_dcd.exists() else traj_path,
            json_path=traj_path.parent / f"{safe}_rep{rep:02d}.json",
            state_csv_path=traj_path.parent / f"{safe}_rep{rep:02d}_md_state.csv",
            mutation=_display_mutation(job["mutation"]),
            replicate=rep,
        )
        job["production_ps"] = float(call.production_ns) * 1000.0

    rows = []
    with ProcessPoolExecutor(max_workers=max(1, args.workers)) as pool:
        futures = {pool.submit(_process_job, j, int(args.pocket_max_frames)): j for j in jobs}
        for fut in as_completed(futures):
            j = futures[fut]
            try:
                rows.append(fut.result())
            except Exception as exc:  # noqa: BLE001
                LOGGER.warning("%s rep%s failed: %s", j["mutation"], j["replicate"], exc)

    per_rep = pd.DataFrame(rows).sort_values(["mutation", "replicate"]).reset_index(drop=True)
    per_rep = per_rep[["mutation", "replicate", "pocket_volume_mean", "pocket_volume_std", "total_ns"]]
    per_rep.to_csv(tables / "pocket_volume_per_rep.csv", index=False)

    per_geno = (
        per_rep.groupby("mutation", as_index=False)["pocket_volume_mean"]
        .agg(volume_mean="mean",
             volume_sem=lambda s: float(s.std(ddof=1) / np.sqrt(len(s))) if len(s) > 1 else 0.0)
    )
    per_geno.to_csv(tables / "pocket_volume_genotype.csv", index=False)

    (cfg / "run_config.json").write_text(json.dumps({
        "manifest": str(args.manifest),
        "pocket_max_frames": int(args.pocket_max_frames),
        "n_replicates": len(per_rep),
    }, indent=2) + "\n")

    LOGGER.info("Wrote %s (%d replicates)", tables / "pocket_volume_per_rep.csv", len(per_rep))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
