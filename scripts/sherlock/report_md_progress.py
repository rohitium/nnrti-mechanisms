#!/usr/bin/env python3
"""Summarize MD extension progress and identify tasks to resume.

This script is intended for Sherlock login-node use.
It reads prepared-system assets + per-replicate JSON outputs and reports:
- how many tasks are complete at target production length
- how many are incomplete / missing
- which tasks are currently running
- which latest logs show segmentation faults
"""

from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
from dataclasses import dataclass
from pathlib import Path


@dataclass
class TaskState:
    mutation: str
    rep: str
    rep_dir: Path
    json_path: Path
    status: str
    steps: int
    running: bool
    latest_log: Path | None
    latest_log_jobid: int | None
    segfault_in_latest_log: bool

    @property
    def key(self) -> str:
        return f"{self.mutation}:rep_{self.rep}"


def _target_steps(ns: float) -> int:
    return max(1, int(round((ns * 1_000_000.0) / 2.0)))


def _safe_steps(payload: dict) -> int:
    value = payload.get("md_production_steps_completed", payload.get("md_production_steps", 0))
    try:
        return int(value or 0)
    except Exception:
        return 0


def _active_md_job_names(user: str) -> set[str]:
    try:
        out = subprocess.check_output(
            ["squeue", "-u", user, "-h", "-t", "PD,R", "-o", "%j"],
            text=True,
            stderr=subprocess.DEVNULL,
        )
    except Exception:
        return set()
    names = set()
    for line in out.splitlines():
        name = line.strip()
        if name.startswith("md_"):
            names.add(name)
    return names


def _latest_logs_by_task(logs_dir: Path) -> dict[tuple[str, str], tuple[int, Path]]:
    latest: dict[tuple[str, str], tuple[int, Path]] = {}
    if not logs_dir.is_dir():
        return latest
    pat = re.compile(r"^md_(.+)_rep([0-9]{2})_([0-9]+)\.log$")
    for path in logs_dir.glob("md_*_rep??_*.log"):
        m = pat.match(path.name)
        if not m:
            continue
        mutation, rep, jobid_s = m.group(1), m.group(2), m.group(3)
        jobid = int(jobid_s)
        key = (mutation, rep)
        prev = latest.get(key)
        if prev is None or jobid > prev[0]:
            latest[key] = (jobid, path)
    return latest


def _has_segfault(path: Path) -> bool:
    try:
        text = path.read_text(errors="ignore")
    except Exception:
        return False
    return "Segmentation fault" in text


def collect_states(root: Path, target_steps: int, user: str) -> list[TaskState]:
    md_runs = root / "results" / "md_runs"
    logs_dir = root / "logs"
    active_names = _active_md_job_names(user)
    latest_logs = _latest_logs_by_task(logs_dir)

    states: list[TaskState] = []
    for system_xml in sorted(md_runs.glob("*/rep_*/assets/*_system.xml")):
        rep_dir = system_xml.parent.parent
        mutation = rep_dir.parent.name
        rep_token = rep_dir.name
        if not rep_token.startswith("rep_"):
            continue
        rep = rep_token.replace("rep_", "", 1)
        json_path = rep_dir / f"{mutation}_rep{rep}.json"

        status = ""
        steps = 0
        if json_path.exists():
            try:
                payload = json.loads(json_path.read_text())
            except Exception:
                payload = {}
            status = str(payload.get("status", "")).lower()
            steps = _safe_steps(payload)

        running = f"md_{mutation}_{rep}" in active_names
        latest = latest_logs.get((mutation, rep))
        latest_log = latest[1] if latest else None
        latest_jobid = latest[0] if latest else None
        segfault = bool(latest_log and _has_segfault(latest_log))

        state = TaskState(
            mutation=mutation,
            rep=rep,
            rep_dir=rep_dir,
            json_path=json_path,
            status=status,
            steps=steps,
            running=running,
            latest_log=latest_log,
            latest_log_jobid=latest_jobid,
            segfault_in_latest_log=segfault,
        )
        states.append(state)

    return states


def main() -> int:
    parser = argparse.ArgumentParser(description="Report MD extension progress")
    parser.add_argument("--root", type=Path, default=Path("."))
    parser.add_argument("--target-ns", type=float, default=float(os.environ.get("MD_PRODUCTION_NS", "10.0")))
    parser.add_argument("--show-incomplete", action="store_true", help="Print per-task lines for incomplete tasks")
    parser.add_argument("--write-lists", action="store_true", help="Write complete/incomplete key lists under results/.status")
    args = parser.parse_args()

    root = args.root.resolve()
    user = os.environ.get("USER", "")
    target_steps = _target_steps(args.target_ns)
    states = collect_states(root, target_steps, user)

    total = len(states)
    complete = sum(1 for s in states if s.status == "ok" and s.steps >= target_steps)
    ok_below = sum(1 for s in states if s.status == "ok" and s.steps < target_steps)
    missing_json = sum(1 for s in states if not s.json_path.exists())
    non_ok_json = sum(1 for s in states if s.json_path.exists() and s.status not in {"", "ok"})
    running = sum(1 for s in states if s.running)
    incomplete_not_running = sum(
        1
        for s in states
        if (s.status != "ok" or s.steps < target_steps) and not s.running
    )
    segfault_latest = sum(
        1
        for s in states
        if (s.status != "ok" or s.steps < target_steps) and s.segfault_in_latest_log
    )

    print("==========================================")
    print("MD Progress Summary")
    print("==========================================")
    print(f"Root:                     {root}")
    print(f"Target ns:                {args.target_ns}")
    print(f"Target steps:             {target_steps}")
    print(f"Prepared systems:         {total}")
    print(f"Complete at target:       {complete}")
    print(f"OK but below target:      {ok_below}")
    print(f"Missing JSON:             {missing_json}")
    print(f"Non-ok JSON:              {non_ok_json}")
    print(f"Currently running (PD/R): {running}")
    print(f"Incomplete not running:   {incomplete_not_running}")
    print(f"Segfault in latest log:   {segfault_latest}")

    if args.write_lists:
        out_dir = root / "results" / ".status"
        out_dir.mkdir(parents=True, exist_ok=True)

        complete_path = out_dir / "complete_at_target.txt"
        incomplete_path = out_dir / "incomplete_not_running.txt"
        with complete_path.open("w") as fh_complete, incomplete_path.open("w") as fh_incomplete:
            for s in sorted(states, key=lambda x: (x.mutation, x.rep)):
                if s.status == "ok" and s.steps >= target_steps:
                    fh_complete.write(f"{s.mutation} {s.rep} {s.rep_dir}\n")
                elif not s.running:
                    fh_incomplete.write(f"{s.mutation} {s.rep} {s.rep_dir}\n")
        print(f"Wrote:                    {complete_path}")
        print(f"Wrote:                    {incomplete_path}")

    if args.show_incomplete:
        print("")
        print("Incomplete tasks:")
        for s in sorted(states, key=lambda x: (x.mutation, x.rep)):
            if s.status == "ok" and s.steps >= target_steps:
                continue
            log_hint = s.latest_log.name if s.latest_log else "-"
            seg = "yes" if s.segfault_in_latest_log else "no"
            run = "yes" if s.running else "no"
            print(
                f"- {s.mutation} rep_{s.rep}: status={s.status or 'missing'} "
                f"steps={s.steps} running={run} segfault_latest={seg} log={log_hint}"
            )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
