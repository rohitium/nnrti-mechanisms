#!/usr/bin/env python3
"""Audit pmx NEQ panel progress from manifest + on-disk artifacts."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import defaultdict
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.fep_pmx.config import FEP_PMX_ROOT


def _neq_root(row: dict[str, str]) -> Path:
    return (
        FEP_PMX_ROOT
        / "legs"
        / row["leg_id"]
        / row["phase"]
        / f"rep_{int(row['replicate']):02d}"
        / "neq"
    )


def _task_ok(row: dict[str, str]) -> tuple[bool, str]:
    neq = _neq_root(row)
    stage = row["stage"]

    if stage == "em":
        path = neq / "em" / "em.gro"
        if path.is_file() and (neq / "em" / "status.json").is_file():
            return True, str(path)
        return False, f"missing {path}"

    if stage == "equil":
        work = neq / row["run_dir"]
        path = work / "equil.gro"
        trr = work / "equil.trr"
        if path.is_file() and trr.is_file() and (work / "status.json").is_file():
            return True, str(trr)
        partial = work / "equil.cpt" if (work / "equil.cpt").is_file() else None
        if partial:
            return False, f"incomplete equil (checkpoint only): {work}"
        if path.is_file() and (work / "status.json").is_file() and not trr.is_file():
            return False, f"missing {trr} (status.json present — re-run equil or delete stale marker)"
        return False, f"missing {path}"

    if stage == "extract":
        snap_dir = neq / row["run_dir"]
        marker = snap_dir / "status.json"
        if not marker.is_file():
            return False, f"missing {marker}"
        prep = json.loads((neq / "neq_prepare.json").read_text())
        n_expected = int(prep["n_snapshots"])
        frames = sorted(snap_dir.glob("frame_*.gro"))
        if len(frames) < n_expected:
            return False, f"{len(frames)}/{n_expected} frames in {snap_dir}"
        return True, f"{len(frames)} frames in {snap_dir}"

    if stage == "switch":
        start = int(row["snapshot_index"])
        end_raw = row.get("snapshot_index_end", "")
        end = int(end_raw) if end_raw not in ("", None) else start
        direction = row["direction"]
        missing: list[str] = []
        for idx in range(start, end + 1):
            dgdl = neq / f"switches/{direction}_{idx:03d}" / "dgdl.xvg"
            if not dgdl.is_file():
                missing.append(str(dgdl))
        if not missing:
            return True, f"{end - start + 1} switches {direction} [{start}-{end}]"
        return False, f"missing {len(missing)}/{end - start + 1} (first: {missing[0]})"

    return False, f"unknown stage {stage}"


def audit_manifest(manifest: Path) -> int:
    rows: list[dict[str, str]] = []
    with manifest.open(newline="") as handle:
        reader = csv.DictReader(handle)
        key = "panel_task_id" if "panel_task_id" in (reader.fieldnames or []) else "task_id"
        for row in reader:
            row = dict(row)
            row["_id"] = row[key]
            rows.append(row)

    by_stage: dict[str, list[tuple[bool, dict[str, str], str]]] = defaultdict(list)
    for row in rows:
        ok, detail = _task_ok(row)
        by_stage[row["stage"]].append((ok, row, detail))

    print(f"Manifest: {manifest}")
    print(f"Total tasks: {len(rows)}")
    print()

    exit_code = 0
    stage_order = ("em", "equil", "extract", "switch")
    for stage in stage_order:
        items = by_stage.get(stage, [])
        if not items:
            continue
        ok_n = sum(1 for ok, _, _ in items if ok)
        print(f"=== {stage.upper()} ({ok_n}/{len(items)} ok) ===")
        if ok_n == len(items):
            print("  all complete")
        else:
            exit_code = 1
            for ok, row, detail in items:
                if ok:
                    continue
                print(
                    f"  FAIL panel={row['_id']} "
                    f"{row['leg_id']} {row['phase']} rep{row['replicate']} "
                    f"→ {detail}"
                )
        print()

    # Per leg-phase-rep rollup: count individual switches (not bundle tasks)
    units: dict[tuple[str, str, int], dict[str, int]] = defaultdict(
        lambda: {"switch_ok": 0, "switch_total": 0}
    )
    for _ok, row, _ in by_stage.get("switch", []):
        key = (row["leg_id"], row["phase"], int(row["replicate"]))
        start = int(row["snapshot_index"])
        end_raw = row.get("snapshot_index_end", "")
        end = int(end_raw) if end_raw not in ("", None) else start
        n = end - start + 1
        units[key]["switch_total"] += n
        direction = row["direction"]
        neq = _neq_root(row)
        for idx in range(start, end + 1):
            if (neq / f"switches/{direction}_{idx:03d}" / "dgdl.xvg").is_file():
                units[key]["switch_ok"] += 1

    if units:
        print("=== SWITCH BY SYSTEM (for BAR) ===")
        for (leg, phase, rep), counts in sorted(units.items()):
            status = "ready" if counts["switch_ok"] == counts["switch_total"] else "incomplete"
            print(
                f"  {leg} {phase} rep{rep:02d}: "
                f"{counts['switch_ok']}/{counts['switch_total']} switches ({status})"
            )
        print()

    return exit_code


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Audit pmx NEQ panel completion.")
    parser.add_argument(
        "--manifest",
        type=Path,
        default=FEP_PMX_ROOT / "neq_panel_manifest.csv",
    )
    args = parser.parse_args(argv)
    if not args.manifest.is_file():
        print(f"Missing manifest: {args.manifest}", file=sys.stderr)
        return 1
    return audit_manifest(args.manifest)


if __name__ == "__main__":
    raise SystemExit(main())
