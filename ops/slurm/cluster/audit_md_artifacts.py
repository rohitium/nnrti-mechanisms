#!/usr/bin/env python3
"""Audit replicate-level MD metadata and compare two repo roots if needed."""
from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from dataclasses import asdict, dataclass
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from nnrti.md.artifact_steps import infer_best_steps, infer_state_csv_path


JSON_PAT = re.compile(r".*_(?:apo_)?rep[0-9]{2}\.json$")


@dataclass(frozen=True)
class ReplicateAudit:
    key: str
    rel_rep_dir: str
    json_path: str
    state_csv_path: str
    status: str
    json_steps: int
    state_csv_steps: int
    best_steps: int
    best_source: str
    json_state_mismatch: bool


def _primary_json(rep_dir: Path) -> Path | None:
    matches = sorted(path for path in rep_dir.glob("*.json") if JSON_PAT.fullmatch(path.name))
    if not matches:
        return None
    return matches[0]


def _status(json_path: Path | None) -> str:
    if json_path is None or not json_path.exists():
        return ""
    try:
        return str(json.loads(json_path.read_text()).get("status", "")).lower()
    except Exception:
        return ""


def _scan_rep_dirs(root: Path, include_apo: bool) -> list[Path]:
    md_runs = root / "results" / "md_runs"
    rep_dirs = sorted(path for path in md_runs.glob("*/rep_*") if path.is_dir())
    if include_apo:
        rep_dirs.extend(sorted(path for path in (md_runs / "apo").glob("*/rep_*") if path.is_dir()))
    return rep_dirs


def audit_root(root: Path, include_apo: bool) -> dict[str, ReplicateAudit]:
    rows: dict[str, ReplicateAudit] = {}
    for rep_dir in _scan_rep_dirs(root, include_apo=include_apo):
        rel_rep_dir = str(rep_dir.relative_to(root))
        json_path = _primary_json(rep_dir)
        state_csv_path = infer_state_csv_path(json_path) if json_path is not None else None
        inferred = infer_best_steps(json_path=json_path, state_csv_path=state_csv_path)
        rows[rel_rep_dir] = ReplicateAudit(
            key=rel_rep_dir,
            rel_rep_dir=rel_rep_dir,
            json_path=str(json_path.relative_to(root)) if json_path is not None and json_path.exists() else "",
            state_csv_path=(
                str(state_csv_path.relative_to(root))
                if state_csv_path is not None and state_csv_path.exists()
                else ""
            ),
            status=_status(json_path),
            json_steps=inferred.json_steps,
            state_csv_steps=inferred.state_csv_steps,
            best_steps=inferred.best_steps,
            best_source=inferred.best_source,
            json_state_mismatch=(
                inferred.json_steps > 0
                and inferred.state_csv_steps > 0
                and inferred.json_steps != inferred.state_csv_steps
            ),
        )
    return rows


def _winner_label(left: ReplicateAudit | None, right: ReplicateAudit | None, left_label: str, right_label: str) -> str:
    if left is None and right is None:
        return "missing"
    if left is None:
        return right_label
    if right is None:
        return left_label
    if left.best_steps > right.best_steps:
        return left_label
    if right.best_steps > left.best_steps:
        return right_label
    if left.status == "ok" and right.status != "ok":
        return left_label
    if right.status == "ok" and left.status != "ok":
        return right_label
    return "tie"


def main() -> int:
    parser = argparse.ArgumentParser(description="Audit MD JSON/state CSV agreement")
    parser.add_argument("--root", type=Path, default=Path("."), help="Primary repo root")
    parser.add_argument("--compare-root", type=Path, default=None, help="Optional second repo root")
    parser.add_argument("--include-apo", action="store_true", help="Include results/md_runs/apo")
    parser.add_argument("--csv-out", type=Path, default=None, help="Optional CSV output path")
    args = parser.parse_args()

    root = args.root.resolve()
    left = audit_root(root, include_apo=args.include_apo)

    right = None
    compare_root = None
    if args.compare_root is not None:
        compare_root = args.compare_root.resolve()
        right = audit_root(compare_root, include_apo=args.include_apo)

    rows = []
    keys = sorted(set(left) | (set(right) if right is not None else set()))
    for key in keys:
        lrow = left.get(key)
        rrow = right.get(key) if right is not None else None
        row = {
            "replicate_key": key,
            "winner": _winner_label(lrow, rrow, "root", "compare_root"),
        }
        if lrow is not None:
            for k, v in asdict(lrow).items():
                row[f"root_{k}"] = v
        if rrow is not None:
            for k, v in asdict(rrow).items():
                row[f"compare_{k}"] = v
        rows.append(row)

    mismatch_count = sum(1 for row in left.values() if row.json_state_mismatch)
    print("Primary root:         {}".format(root))
    print("Replicates audited:   {}".format(len(left)))
    print("JSON/state mismatch:  {}".format(mismatch_count))
    if compare_root is not None and right is not None:
        print("Compare root:         {}".format(compare_root))
        print("Compare audited:      {}".format(len(right)))
        print("Root wins:            {}".format(sum(1 for row in rows if row["winner"] == "root")))
        print("Compare wins:         {}".format(sum(1 for row in rows if row["winner"] == "compare_root")))
        print("Ties:                 {}".format(sum(1 for row in rows if row["winner"] == "tie")))

    if args.csv_out is not None:
        args.csv_out.parent.mkdir(parents=True, exist_ok=True)
        fieldnames = sorted({k for row in rows for k in row})
        with args.csv_out.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows)
        print("Wrote CSV:            {}".format(args.csv_out))
    else:
        for row in rows:
            if right is None:
                if row.get("root_json_state_mismatch"):
                    print(
                        "- {}: status={} best_steps={} source={} json_steps={} state_csv_steps={}".format(
                            row["replicate_key"],
                            row.get("root_status", ""),
                            row.get("root_best_steps", 0),
                            row.get("root_best_source", ""),
                            row.get("root_json_steps", 0),
                            row.get("root_state_csv_steps", 0),
                        )
                    )
            else:
                if row["winner"] != "tie":
                    print(
                        "- {}: winner={} root_steps={} compare_steps={}".format(
                            row["replicate_key"],
                            row["winner"],
                            row.get("root_best_steps", 0),
                            row.get("compare_best_steps", 0),
                        )
                    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
