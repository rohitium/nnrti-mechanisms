#!/usr/bin/env python3
"""Run apo OpenMM MD for genotypes with prepared assets but no trajectories."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[3]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from nnrti.fep.mutations import apo_safe_label, canonical_label, safe_label


def _apo_tasks(
    mutations: list[str],
    *,
    apo_root: Path,
    holo_root: Path,
    replicates: tuple[int, ...],
) -> list[dict]:
    tasks: list[dict] = []
    for mutation in mutations:
        label = canonical_label(mutation)
        apo_dir = apo_safe_label(label)
        holo_dir = safe_label(label)
        for rep in replicates:
            apo_rep = apo_root / apo_dir / f"rep_{rep:02d}"
            assets = apo_rep / "assets"
            topology = assets / f"{apo_dir}_apo_md_rep{rep:02d}_start.pdb"
            system_xml = assets / f"{apo_dir}_apo_md_rep{rep:02d}_system.xml"
            output_json = apo_rep / f"{apo_dir}_apo_rep{rep:02d}.json"
            minimized = holo_root / holo_dir / f"rep_{rep:02d}" / f"{holo_dir}_minimized_rep{rep:02d}.pdb"
            tasks.append(
                {
                    "mutation": label,
                    "safe_label": apo_dir,
                    "replicate": rep,
                    "topology_pdb": topology,
                    "system_xml": system_xml,
                    "output_json": output_json,
                    "minimized_pdb": minimized,
                }
            )
    return tasks


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Run apo MD for prepared genotypes.")
    parser.add_argument(
        "--mutations",
        nargs="+",
        default=["Y188L"],
        help="Mutation labels (default: Y188L P0 blocker).",
    )
    parser.add_argument("--replicates", type=int, nargs="+", default=[1, 2, 3])
    parser.add_argument("--apo-root", type=Path, default=Path("results/md_runs/apo"))
    parser.add_argument("--holo-root", type=Path, default=Path("results/md_runs"))
    parser.add_argument("--production-ns", type=float, default=100.0)
    parser.add_argument("--heating-ps", type=float, default=25.0)
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args(argv)

    from nnrti.md.manifest import MDTask
    from nnrti.md.worker import run_md_task

    tasks = _apo_tasks(
        args.mutations,
        apo_root=args.apo_root,
        holo_root=args.holo_root,
        replicates=tuple(args.replicates),
    )
    pending = []
    for spec in tasks:
        missing = [p for p in (spec["topology_pdb"], spec["system_xml"]) if not p.is_file()]
        if missing:
            print(f"SKIP missing: {missing[0]}", file=sys.stderr)
            continue
        pending.append(spec)

    if args.dry_run:
        print(json.dumps([{k: str(v) for k, v in t.items()} for t in pending], indent=2))
        return 0

    if not pending:
        print("No runnable apo tasks.", file=sys.stderr)
        return 1

    for idx, spec in enumerate(pending):
        task = MDTask(
            task_id=idx,
            structure="DOR",
            mutation=spec["mutation"],
            safe_label=spec["safe_label"],
            replicate=spec["replicate"],
            minimized_pdb=str(spec["minimized_pdb"]),
            ligand_sdf="",
            ligand_resname="",
            fold_reduction=None,
            output_json=str(spec["output_json"]),
            leg="apo",
            prepared_topology_pdb=str(spec["topology_pdb"]),
            prepared_system_xml=str(spec["system_xml"]),
        )
        print(f"Running apo MD: {spec['mutation']} rep {spec['replicate']:02d}")
        result = run_md_task(
            task,
            heating_ps=args.heating_ps,
            production_ns=args.production_ns,
            report_interval=2000,
            checkpoint_interval=5000,
            resume_from_checkpoint=True,
            force=args.force,
        )
        print(f"  status={result.get('status')}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
