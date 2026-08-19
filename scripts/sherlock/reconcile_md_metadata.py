#!/usr/bin/env python3
"""Reconcile MD JSON metadata with paired md_state.csv files."""
from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.md.artifact_steps import infer_state_csv_path, reconcile_json_with_state_csv


JSON_PAT = re.compile(r".*_(?:apo_)?rep[0-9]{2}\.json$")


def _scan_json_paths(root: Path, include_apo: bool) -> list[Path]:
    md_runs = root / "results" / "md_runs"
    paths = sorted(path for path in md_runs.glob("*/*/*.json") if JSON_PAT.fullmatch(path.name))
    if include_apo:
        paths.extend(sorted(path for path in (md_runs / "apo").glob("*/*/*.json") if JSON_PAT.fullmatch(path.name)))
    return paths


def _target_steps_from_ns(ns: float | None) -> int | None:
    if ns is None:
        return None
    return max(1, int(round((float(ns) * 1_000_000.0) / 2.0)))


def main() -> int:
    parser = argparse.ArgumentParser(description="Reconcile MD JSON files with md_state.csv")
    parser.add_argument("--root", type=Path, default=Path("."))
    parser.add_argument("--include-apo", action="store_true")
    parser.add_argument("--target-ns", type=float, default=None)
    parser.add_argument("--write", action="store_true", help="Rewrite stale JSON files in place")
    parser.add_argument(
        "--allow-downgrade",
        action="store_true",
        help=(
            "Permit LOWERING md_production_steps_completed from state.csv. Off by "
            "default: state.csv is often a stale mid-slice dump, and a downgrade "
            "silently compresses every analysis time axis. Prefer the analysis-DCD "
            "fingerprint (src.analysis.md_timing) before using this."
        ),
    )
    parser.add_argument("--csv-out", type=Path, default=None)
    args = parser.parse_args()

    root = args.root.resolve()
    target_steps = _target_steps_from_ns(args.target_ns)
    rows: list[dict] = []
    changed = 0
    mismatched = 0

    for json_path in _scan_json_paths(root, include_apo=args.include_apo):
        state_csv_path = infer_state_csv_path(json_path)
        before = reconcile_json_with_state_csv(
            json_path,
            state_csv_path=state_csv_path,
            write=False,
            target_steps=target_steps,
        )
        after = reconcile_json_with_state_csv(
            json_path,
            state_csv_path=state_csv_path,
            write=args.write,
            allow_downgrade=args.allow_downgrade,
            target_steps=target_steps,
        )
        if not before.consistent:
            mismatched += 1
        if after.changed:
            changed += 1
        rows.append(
            {
                "json_path": str(json_path.relative_to(root)),
                "state_csv_path": (
                    str(state_csv_path.relative_to(root))
                    if state_csv_path is not None and state_csv_path.exists()
                    else ""
                ),
                "status": after.status,
                "json_steps_before": before.json_steps,
                "state_csv_steps": before.state_csv_steps,
                "consistent_before": before.consistent,
                "json_steps_after": after.json_steps,
                "consistent_after": after.consistent,
                "changed": after.changed,
            }
        )

    print("Root:                 {}".format(root))
    print("JSON files scanned:   {}".format(len(rows)))
    print("Mismatched before:    {}".format(mismatched))
    print("JSON files changed:   {}".format(changed))
    print("Write mode:           {}".format("yes" if args.write else "no"))
    if target_steps is not None:
        print("Target steps:         {}".format(target_steps))

    if args.csv_out is not None:
        args.csv_out.parent.mkdir(parents=True, exist_ok=True)
        with args.csv_out.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()) if rows else [])
            if rows:
                writer.writeheader()
                writer.writerows(rows)
        print("Wrote CSV:            {}".format(args.csv_out))
    else:
        for row in rows:
            if not row["consistent_before"] or row["changed"]:
                print(
                    "- {}: json_before={} state_csv={} json_after={} changed={}".format(
                        row["json_path"],
                        row["json_steps_before"],
                        row["state_csv_steps"],
                        row["json_steps_after"],
                        row["changed"],
                    )
                )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
