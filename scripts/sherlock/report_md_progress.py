#!/usr/bin/env python3
"""Summarize MD extension progress and identify tasks to resume."""
from __future__ import annotations

import argparse
import glob
import json
import os
import re
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from nnrti.md.artifact_steps import infer_state_csv_path, reconcile_json_with_state_csv


class TaskState(object):
    def __init__(
        self,
        mutation,
        rep,
        rep_dir,
        json_path,
        status,
        steps,
        json_steps,
        state_csv_steps,
        steps_consistent,
        running,
        latest_log,
        latest_log_jobid,
        segfault_in_latest_log,
    ):
        self.mutation = mutation
        self.rep = rep
        self.rep_dir = rep_dir
        self.json_path = json_path
        self.status = status
        self.steps = steps
        self.json_steps = json_steps
        self.state_csv_steps = state_csv_steps
        self.steps_consistent = steps_consistent
        self.running = running
        self.latest_log = latest_log
        self.latest_log_jobid = latest_log_jobid
        self.segfault_in_latest_log = segfault_in_latest_log


def _target_steps(ns):
    return max(1, int(round((float(ns) * 1_000_000.0) / 2.0)))


def _active_md_job_names(user):
    try:
        out = subprocess.check_output(
            ["squeue", "-u", user, "-h", "-t", "PD,R", "-o", "%j"],
            universal_newlines=True,
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


def _latest_logs_by_task(logs_dir):
    latest = {}
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


def _has_segfault(path):
    try:
        text = path.read_text(errors="ignore")
    except Exception:
        return False
    return "Segmentation fault" in text


def collect_states(root, user, repair=False, target_steps=None):
    md_runs = root / "results" / "md_runs"
    logs_dir = root / "logs"
    active_names = _active_md_job_names(user)
    latest_logs = _latest_logs_by_task(logs_dir)

    states = []
    for system_xml in sorted(glob.glob(str(md_runs / "*/rep_*/assets/*_system.xml"))):
        system_xml = Path(system_xml)
        rep_dir = system_xml.parent.parent
        mutation = rep_dir.parent.name
        rep_token = rep_dir.name
        if not rep_token.startswith("rep_"):
            continue
        rep = rep_token.replace("rep_", "", 1)
        json_path = rep_dir / "{}_rep{}.json".format(mutation, rep)

        reconciled = reconcile_json_with_state_csv(
            json_path=json_path,
            state_csv_path=infer_state_csv_path(json_path),
            write=bool(repair),
            target_steps=target_steps,
        )

        running = "md_{}_{}".format(mutation, rep) in active_names
        latest = latest_logs.get((mutation, rep))
        latest_log = latest[1] if latest else None
        latest_jobid = latest[0] if latest else None
        segfault = bool(latest_log and _has_segfault(latest_log))

        states.append(
            TaskState(
                mutation=mutation,
                rep=rep,
                rep_dir=rep_dir,
                json_path=json_path,
                status=reconciled.status,
                steps=reconciled.json_steps,
                json_steps=reconciled.json_steps,
                state_csv_steps=reconciled.state_csv_steps,
                steps_consistent=reconciled.consistent,
                running=running,
                latest_log=latest_log,
                latest_log_jobid=latest_jobid,
                segfault_in_latest_log=segfault,
            )
        )

    return states


def main():
    parser = argparse.ArgumentParser(description="Report MD extension progress")
    parser.add_argument("--root", type=Path, default=Path("."))
    parser.add_argument(
        "--target-ns",
        type=float,
        default=float(os.environ.get("MD_PRODUCTION_NS", "100.0")),
    )
    parser.add_argument(
        "--show-incomplete",
        action="store_true",
        help="Print per-task lines for incomplete tasks",
    )
    parser.add_argument(
        "--write-lists",
        action="store_true",
        help="Write complete/incomplete key lists under results/.status",
    )
    parser.add_argument(
        "--repair",
        action="store_true",
        help="Rewrite stale JSON step metadata from md_state.csv before reporting",
    )
    args = parser.parse_args()

    root = args.root.resolve()
    user = os.environ.get("USER", "")
    target_steps = _target_steps(args.target_ns)
    states = collect_states(root, user, repair=args.repair, target_steps=target_steps)

    total = len(states)
    complete = sum(1 for s in states if s.status == "ok" and s.steps >= target_steps)
    ok_below = sum(1 for s in states if s.status == "ok" and s.steps < target_steps)
    missing_json = sum(1 for s in states if not s.json_path.exists())
    non_ok_json = sum(1 for s in states if s.json_path.exists() and s.status not in {"", "ok"})
    mismatched = sum(1 for s in states if not s.steps_consistent)
    running = sum(1 for s in states if s.running)
    incomplete_not_running = sum(
        1 for s in states if (s.status != "ok" or s.steps < target_steps) and not s.running
    )
    segfault_latest = sum(
        1
        for s in states
        if (s.status != "ok" or s.steps < target_steps) and s.segfault_in_latest_log
    )

    print("==========================================")
    print("MD Progress Summary")
    print("==========================================")
    print("Root:                     {}".format(root))
    print("Target ns:                {}".format(args.target_ns))
    print("Target steps:             {}".format(target_steps))
    print("Prepared systems:         {}".format(total))
    print("Complete at target:       {}".format(complete))
    print("OK but below target:      {}".format(ok_below))
    print("Missing JSON:             {}".format(missing_json))
    print("Non-ok JSON:              {}".format(non_ok_json))
    print("JSON/state mismatches:    {}".format(mismatched))
    print("Currently running (PD/R): {}".format(running))
    print("Incomplete not running:   {}".format(incomplete_not_running))
    print("Segfault in latest log:   {}".format(segfault_latest))

    if args.write_lists:
        out_dir = root / "results" / ".status"
        out_dir.mkdir(parents=True, exist_ok=True)
        complete_path = out_dir / "complete_at_target.txt"
        incomplete_path = out_dir / "incomplete_not_running.txt"

        with complete_path.open("w") as fh_complete, incomplete_path.open("w") as fh_incomplete:
            for s in sorted(states, key=lambda x: (x.mutation, x.rep)):
                if s.status == "ok" and s.steps >= target_steps:
                    fh_complete.write("{} {} {}\n".format(s.mutation, s.rep, s.rep_dir))
                elif not s.running:
                    fh_incomplete.write("{} {} {}\n".format(s.mutation, s.rep, s.rep_dir))

        print("Wrote:                    {}".format(complete_path))
        print("Wrote:                    {}".format(incomplete_path))

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
                "- {} rep_{}: status={} steps={} consistent={} json_steps={} state_csv_steps={} running={} segfault_latest={} log={}".format(
                    s.mutation,
                    s.rep,
                    s.status or "missing",
                    s.steps,
                    "yes" if s.steps_consistent else "no",
                    s.json_steps,
                    s.state_csv_steps,
                    run,
                    seg,
                    log_hint,
                )
            )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
