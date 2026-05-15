#!/usr/bin/env python3
"""Compare MM/GBSA estimates across trajectory sampling windows."""
from __future__ import annotations

import argparse
import concurrent.futures as cf
from pathlib import Path

import MDAnalysis as mda
import pandas as pd

from src.analysis.cli.compute_mmgbsa_safe import (
    _infer_rep_dir,
    _infer_total_steps,
    _nonempty_path,
    _resolve_local_path,
)
from src.analysis.result_collector import collect_md_results
from src.md.openmm.mmgbsa import _select_snapshot_indices, compute_mmgbsa_from_trajectory


def _task_from_payload(payload: tuple[dict, float, int, str]) -> dict:
    return _task(*payload)


def _task(row_dict: dict, window_ns: float, snapshots: int, ligand_resname: str) -> dict:
    row = pd.Series(row_dict)
    mutation = str(row["mutation"])
    replicate = int(row["replicate"])
    safe = str(row["safe_label"])
    rep_dir = _infer_rep_dir(row)

    min_pdb = _resolve_local_path(
        _nonempty_path(row.get("minimized_pdb")),
        rep_dir / f"{safe}_minimized_rep{replicate:02d}.pdb",
    )
    dcd = _resolve_local_path(
        _nonempty_path(row.get("analysis_dcd")),
        rep_dir / f"{safe}_rep{replicate:02d}_analysis.dcd",
    )
    analysis_topo = _resolve_local_path(
        _nonempty_path(row.get("analysis_topology_pdb")),
        rep_dir / f"{safe}_rep{replicate:02d}_analysis_topology.pdb",
    )
    ligand_sdf = _resolve_local_path(_nonempty_path(row.get("ligand_sdf")))
    if min_pdb is None or dcd is None or analysis_topo is None or ligand_sdf is None:
        raise FileNotFoundError(f"Missing inputs for {mutation} rep{replicate}")

    total_steps = _infer_total_steps(row, rep_dir, safe, replicate)
    total_time_ns = float(total_steps) * 2.0 / 1_000_000.0 if total_steps else None

    universe = mda.Universe(str(analysis_topo), str(dcd))
    n_frames = len(universe.trajectory)
    selected = _select_snapshot_indices(
        n_frames=n_frames,
        discard_fraction=0.25,
        n_snapshots=snapshots,
        dt_ps=getattr(universe.trajectory, "dt", None),
        sample_window_ns=window_ns,
        total_time_ns=total_time_ns,
    )

    result = compute_mmgbsa_from_trajectory(
        minimized_pdb_path=min_pdb,
        trajectory_dcd_path=dcd,
        ligand_resname=ligand_resname,
        ligand_sdf=ligand_sdf,
        n_snapshots=snapshots,
        discard_fraction=0.25,
        sample_window_ns=window_ns,
        total_time_ns=total_time_ns,
        analysis_topology_pdb_path=analysis_topo,
    )

    return {
        "mutation": mutation,
        "safe_label": safe,
        "replicate": replicate,
        "window_ns": float(window_ns),
        "n_frames_total": int(n_frames),
        "total_time_ns": total_time_ns,
        "frame_spacing_ns": (
            float(total_time_ns) / float(n_frames - 1)
            if total_time_ns is not None and n_frames > 1
            else float("nan")
        ),
        "selected_frames": int(len(selected)),
        "first_selected_frame": int(selected[0]),
        "last_selected_frame": int(selected[-1]),
        "binding_dg": result.binding_dg_mean,
        "binding_dg_std": result.binding_dg_std,
        "binding_dg_vdw": result.delta_e_vdw_mean,
        "binding_dg_electrostatic": result.delta_e_elec_mean,
        "binding_dg_gb": result.delta_g_gb_mean,
        "binding_dg_sa": result.delta_g_sa_mean,
        "mmgbsa_snapshots": result.n_snapshots,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, default=Path("manifests/md_manifest.csv"))
    parser.add_argument("--results-dir", type=Path, default=Path("results"))
    parser.add_argument("--mutation", type=str, default="WT", help="Mutation label or comma-separated labels.")
    parser.add_argument("--windows", type=str, default="1,5,10,20")
    parser.add_argument("--snapshots", type=int, default=100)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--ligand-resname", type=str, default="2KW")
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("results/analysis/binding_energy/last1ns/wt_window_sensitivity.csv"),
    )
    args = parser.parse_args()

    md_df = collect_md_results(args.manifest, args.results_dir)
    mutations = {x.strip() for x in str(args.mutation).split(",") if x.strip()}
    target = md_df[md_df["mutation"].astype(str).isin(mutations)].copy()
    if target.empty:
        raise ValueError(f"No rows found for mutation(s) {sorted(mutations)!r}")

    windows = [float(x) for x in str(args.windows).split(",") if x.strip()]
    jobs = [
        (row.to_dict(), window, int(args.snapshots), str(args.ligand_resname))
        for _, row in target.iterrows()
        for window in windows
    ]

    rows: list[dict] = []
    workers = max(1, int(args.workers))
    if workers == 1:
        for job in jobs:
            rows.append(_task(*job))
    else:
        with cf.ProcessPoolExecutor(max_workers=workers) as pool:
            for row in pool.map(_task_from_payload, jobs):
                rows.append(row)

    out = pd.DataFrame(rows).sort_values(["replicate", "window_ns"], kind="stable")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.output, index=False)
    print(out.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
