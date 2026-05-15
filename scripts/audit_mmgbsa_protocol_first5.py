from __future__ import annotations

import argparse
import concurrent.futures as cf
import os
from pathlib import Path

import pandas as pd

from src.analysis.cli.compute_mmgbsa_safe import (
    _compute_one_task,
    _infer_rep_dir,
    _nonempty_path,
    _resolve_local_path,
)
from src.analysis.result_collector import collect_md_results


ENERGY_COLUMNS = [
    "binding_dg",
    "binding_dg_vdw",
    "binding_dg_electrostatic",
    "binding_dg_gb",
    "binding_dg_sa",
]


def _build_first_tasks(args: argparse.Namespace) -> list[dict]:
    md_df = collect_md_results(args.manifest, args.results_dir)
    if md_df.empty:
        raise RuntimeError("No MD results found.")

    tasks: list[dict] = []
    for _, row in md_df.iterrows():
        rep_dir = _infer_rep_dir(row)
        safe = str(row["safe_label"])
        rep = int(row["replicate"])

        min_pdb = _resolve_local_path(
            _nonempty_path(row.get("minimized_pdb")),
            rep_dir / f"{safe}_minimized_rep{rep:02d}.pdb",
        )
        dcd = _resolve_local_path(
            _nonempty_path(row.get("analysis_dcd")),
            rep_dir / f"{safe}_rep{rep:02d}_analysis.dcd",
        )
        analysis_topo = _resolve_local_path(
            _nonempty_path(row.get("analysis_topology_pdb")),
            rep_dir / f"{safe}_rep{rep:02d}_analysis_topology.pdb",
        )
        ligand_sdf = _resolve_local_path(_nonempty_path(row.get("ligand_sdf")))
        if None in {min_pdb, dcd, analysis_topo, ligand_sdf}:
            continue
        if not min_pdb.exists() or not dcd.exists() or not analysis_topo.exists() or not ligand_sdf.exists():
            continue

        tasks.append(
            {
                "structure": row["structure"],
                "mutation": str(row["mutation"]),
                "safe_label": safe,
                "replicate": rep,
                "fold_reduction": row.get("fold_reduction"),
                "min_pdb": str(min_pdb),
                "dcd": str(dcd),
                "analysis_topo": str(analysis_topo),
                "ligand_sdf": str(ligand_sdf),
                "ligand_resname": args.ligand_resname,
                "snapshots": int(args.snapshots),
                "discard_fraction": float(args.discard_fraction),
                "sample_window_ns": None,
                "total_time_ns": None,
                "time_source": "trajectory_dt",
                "sample_last_frames": None,
            }
        )
        if len(tasks) >= int(args.n):
            break
    return tasks


def _run_tasks(tasks: list[dict], workers: int) -> pd.DataFrame:
    rows: list[dict] = []
    if workers == 1:
        for task in tasks:
            ok, row, err = _compute_one_task(task)
            if not ok or row is None:
                raise RuntimeError(err)
            rows.append(row)
    else:
        with cf.ProcessPoolExecutor(max_workers=workers) as ex:
            futures = {ex.submit(_compute_one_task, task): task for task in tasks}
            for fut in cf.as_completed(futures):
                ok, row, err = fut.result()
                if not ok or row is None:
                    raise RuntimeError(err)
                rows.append(row)
    return pd.DataFrame(rows).sort_values(["mutation", "replicate"], kind="stable").reset_index(drop=True)


def _run_tasks_flush(tasks: list[dict], output_path: Path) -> pd.DataFrame:
    rows: list[dict] = []
    completed: set[tuple[str, int]] = set()
    if output_path.exists():
        existing = pd.read_csv(output_path)
        if not existing.empty:
            rows.extend(existing.to_dict("records"))
            completed = set(zip(existing["mutation"].astype(str), existing["replicate"].astype(int)))

    for i, task in enumerate(tasks, 1):
        key = (str(task["mutation"]), int(task["replicate"]))
        if key in completed:
            print(f"[{i}/{len(tasks)}] skipping {key[0]} rep{key[1]} (already written)", flush=True)
            continue
        print(f"[{i}/{len(tasks)}] computing {key[0]} rep{key[1]}", flush=True)
        ok, row, err = _compute_one_task(task)
        if not ok or row is None:
            raise RuntimeError(err)
        rows.append(row)
        pd.DataFrame(rows).sort_values(["mutation", "replicate"], kind="stable").to_csv(output_path, index=False)
        print(
            f"[{i}/{len(tasks)}] wrote {key[0]} rep{key[1]} "
            f"DG={row['binding_dg']:.3f} vdw={row['binding_dg_vdw']:.3f} gb={row['binding_dg_gb']:.3f}",
            flush=True,
        )

    return pd.DataFrame(rows).sort_values(["mutation", "replicate"], kind="stable").reset_index(drop=True)


def _load_reference(path: Path, label: str) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    df = pd.read_csv(path)
    keep = ["mutation", "replicate", *[c for c in ENERGY_COLUMNS if c in df.columns]]
    out = df[keep].copy()
    out = out.rename(columns={c: f"{label}_{c}" for c in ENERGY_COLUMNS if c in out.columns})
    return out


def _compare(audit: pd.DataFrame, refs: list[tuple[Path, str]]) -> pd.DataFrame:
    cols = ["mutation", "replicate", *ENERGY_COLUMNS, "mmgbsa_snapshots"]
    merged = audit[[c for c in cols if c in audit.columns]].copy()
    for path, label in refs:
        ref = _load_reference(path, label)
        if ref.empty:
            continue
        merged = merged.merge(ref, on=["mutation", "replicate"], how="left")
        for col in ENERGY_COLUMNS:
            ref_col = f"{label}_{col}"
            if ref_col in merged.columns and col in merged.columns:
                merged[f"{label}_minus_audit_{col}"] = merged[ref_col] - merged[col]
    return merged


def main() -> int:
    parser = argparse.ArgumentParser(description="Audit first N MM/GBSA trajectories using the original protocol.")
    parser.add_argument("--manifest", type=Path, default=Path("manifests/md_manifest.csv"))
    parser.add_argument("--results-dir", type=Path, default=Path("results"))
    parser.add_argument("--output-dir", type=Path, default=Path("results/analysis/binding_energy/audit_first5_protocol"))
    parser.add_argument("--n", type=int, default=5)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--snapshots", type=int, default=100)
    parser.add_argument("--discard-fraction", type=float, default=0.25)
    parser.add_argument("--ligand-resname", type=str, default="2KW")
    args = parser.parse_args()

    os.environ.setdefault("OPENMM_PLATFORM", "CPU")
    os.environ.setdefault("OPENMM_CPU_THREADS", "1")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    tasks = _build_first_tasks(args)
    if len(tasks) < args.n:
        raise RuntimeError(f"Only found {len(tasks)} runnable tasks; requested {args.n}.")

    run_label = f"workers{max(1, int(args.workers))}"
    audit_path = args.output_dir / f"mmgbsa_first{args.n}_{run_label}.csv"
    if max(1, int(args.workers)) == 1:
        audit = _run_tasks_flush(tasks, audit_path)
    else:
        audit = _run_tasks(tasks, max(1, int(args.workers)))
        audit.to_csv(audit_path, index=False)

    comparison = _compare(
        audit,
        [
            (Path("results/archive/2026-05-13_binding_energy_pre_recompute/mmgbsa_replicate_metrics.csv"), "old_archive"),
            (Path("results/archive/2026-05-13_binding_energy_pre_recompute/checkpoint_mmgbsa_replicate_metrics.csv"), "archived_checkpoint"),
            (Path("results/mmgbsa_replicate_metrics.csv"), "current_top"),
        ],
    )
    comparison_path = args.output_dir / f"comparison_first{args.n}_{run_label}.csv"
    comparison.to_csv(comparison_path, index=False)

    print(f"Wrote {audit_path}")
    print(f"Wrote {comparison_path}")
    print(comparison[["mutation", "replicate", "binding_dg", "binding_dg_vdw", "binding_dg_gb"]].to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
