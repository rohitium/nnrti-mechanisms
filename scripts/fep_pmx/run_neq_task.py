#!/usr/bin/env python3
"""Run one NEQ manifest task (EM, equil, snapshot extract, or switch)."""

from __future__ import annotations

import argparse
import csv
import json
import os
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.fep_pmx.gromacs_utils import GromacsError, find_gmx, gromacs_env, run_gmx


def _neq_root(row: dict[str, str]) -> Path:
    return (
        REPO_ROOT
        / "results/analysis/fep_pmx/legs"
        / row["leg_id"]
        / row["phase"]
        / f"rep_{int(row['replicate']):02d}"
        / "neq"
    )


def _load_manifest_row(manifest: Path, task_id: int) -> dict[str, str]:
    with manifest.open(newline="") as handle:
        for row in csv.DictReader(handle):
            key = "panel_task_id" if "panel_task_id" in row else "task_id"
            if int(row[key]) == task_id:
                return row
    raise KeyError(f"task_id {task_id} not found in {manifest}")


def _mdrun(gmx_mdrun: str, *, cwd: Path, deffnm: str, env: dict[str, str]) -> None:
    run_gmx(
        gmx_mdrun,
        ["mdrun", "-v", "-deffnm", deffnm, "-nb", "gpu", "-gpu_id", "0"],
        cwd=cwd,
        env=env,
    )


def _run_em(neq: Path, env: dict[str, str], gmx: str, gmx_mdrun: str) -> None:
    work = neq / "em"
    work.mkdir(parents=True, exist_ok=True)
    if (work / "em.gro").is_file() and (work / "status.json").is_file():
        return
    run_gmx(
        gmx,
        [
            "grompp",
            "-f",
            "../mdp/em.mdp",
            "-c",
            "../system.gro",
            "-p",
            "../system.top",
            "-o",
            "em.tpr",
            "-maxwarn",
            "10",
        ],
        cwd=work,
        env=env,
    )
    _mdrun(gmx_mdrun, cwd=work, deffnm="em", env=env)
    (work / "status.json").write_text(json.dumps({"stage": "em", "status": "ok"}) + "\n")


def _run_equil(neq: Path, row: dict[str, str], env: dict[str, str], gmx: str, gmx_mdrun: str) -> None:
    lambda_state = int(row["lambda_state"])
    work = neq / row["run_dir"]
    work.mkdir(parents=True, exist_ok=True)
    if (work / "equil.gro").is_file() and (work / "status.json").is_file():
        return

    em_gro = neq / "em" / "em.gro"
    if not em_gro.is_file():
        raise FileNotFoundError(f"Missing minimized structure: {em_gro}")

    mdp = neq / "mdp" / f"npt_eq_lambda{lambda_state}.mdp"
    run_gmx(
        gmx,
        [
            "grompp",
            "-f",
            f"../mdp/npt_eq_lambda{lambda_state}.mdp",
            "-c",
            "../em/em.gro",
            "-p",
            "../system.top",
            "-o",
            "equil.tpr",
            "-maxwarn",
            "10",
        ],
        cwd=work,
        env=env,
    )
    _mdrun(gmx_mdrun, cwd=work, deffnm="equil", env=env)
    (work / "status.json").write_text(
        json.dumps({"stage": "equil", "lambda_state": lambda_state, "status": "ok"}) + "\n"
    )


def _run_extract(neq: Path, row: dict[str, str], env: dict[str, str], gmx: str) -> None:
    lambda_state = int(row["lambda_state"])
    eq_dir = neq / f"eq_lambda{lambda_state}"
    snap_dir = neq / row["run_dir"]
    snap_dir.mkdir(parents=True, exist_ok=True)
    marker = snap_dir / "status.json"
    if marker.is_file():
        return

    trr = eq_dir / "equil.trr"
    tpr = eq_dir / "equil.tpr"
    if not trr.is_file() or not tpr.is_file():
        raise FileNotFoundError(f"Missing equil trajectory in {eq_dir}")

    meta = json.loads((neq / "neq_prepare.json").read_text())
    for idx, time_ps in enumerate(meta["snapshot_times_ps"]):
        out_gro = snap_dir / f"frame_{idx:03d}.gro"
        if out_gro.is_file():
            continue
        run_gmx(
            gmx,
            [
                "trjconv",
                "-f",
                f"../eq_lambda{lambda_state}/equil.trr",
                "-s",
                f"../eq_lambda{lambda_state}/equil.tpr",
                "-dump",
                str(time_ps),
                "-o",
                out_gro.name,
                "-pbc",
                "mol",
                "-ur",
                "compact",
            ],
            cwd=snap_dir,
            input_text="System\n",
            env=env,
        )
        sidecar = snap_dir / f"frame_{idx:03d}.json"
        sidecar.write_text(json.dumps({"time_ps": time_ps, "lambda_state": lambda_state}) + "\n")

    marker.write_text(json.dumps({"stage": "extract", "lambda_state": lambda_state, "status": "ok"}) + "\n")


def _run_switch(neq: Path, row: dict[str, str], env: dict[str, str], gmx: str, gmx_mdrun: str) -> None:
    direction = row["direction"]
    snapshot_index = int(row["snapshot_index"])
    lambda_state = int(row["lambda_state"])
    work = neq / row["run_dir"]
    work.mkdir(parents=True, exist_ok=True)
    if (work / "dgdl.xvg").is_file() and (work / "status.json").is_file():
        return

    frame_gro = neq / f"snapshots/lambda{lambda_state}" / f"frame_{snapshot_index:03d}.gro"
    if not frame_gro.is_file():
        raise FileNotFoundError(f"Missing snapshot: {frame_gro}")

    frame_rel = Path("..") / "snapshots" / f"lambda{lambda_state}" / f"frame_{snapshot_index:03d}.gro"
    trr_rel = Path("..") / f"eq_lambda{lambda_state}" / "equil.trr"
    mdp_name = "nonequil_fwd.mdp" if direction == "fwd" else "nonequil_rev.mdp"
    sidecar = neq / "snapshots" / f"lambda{lambda_state}" / f"frame_{snapshot_index:03d}.json"
    time_ps = json.loads(sidecar.read_text())["time_ps"]

    run_gmx(
        gmx,
        [
            "grompp",
            "-f",
            f"../mdp/{mdp_name}",
            "-c",
            frame_rel.as_posix(),
            "-t",
            trr_rel.as_posix(),
            "-time",
            str(time_ps),
            "-p",
            "../system.top",
            "-o",
            "switch.tpr",
            "-maxwarn",
            "10",
        ],
        cwd=work,
        env=env,
    )
    _mdrun(gmx_mdrun, cwd=work, deffnm="switch", env=env)
    dhdl = work / "switch.dhdl.xvg"
    if dhdl.is_file() and not (work / "dgdl.xvg").exists():
        dhdl.rename(work / "dgdl.xvg")
    (work / "status.json").write_text(
        json.dumps(
            {
                "stage": "switch",
                "direction": direction,
                "snapshot_index": snapshot_index,
                "time_ps": time_ps,
                "status": "ok",
            }
        )
        + "\n"
    )


def run_task(manifest: Path, task_id: int) -> None:
    row = _load_manifest_row(manifest, task_id)
    neq = _neq_root(row)
    env = gromacs_env()
    gmx = find_gmx()
    gmx_mdrun = env.get("GMX_MDRUN") or gmx

    stage = row["stage"]
    if stage == "em":
        _run_em(neq, env, gmx, gmx_mdrun)
    elif stage == "equil":
        _run_equil(neq, row, env, gmx, gmx_mdrun)
    elif stage == "extract":
        _run_extract(neq, row, env, gmx)
    elif stage == "switch":
        _run_switch(neq, row, env, gmx, gmx_mdrun)
    else:
        raise ValueError(f"Unknown stage: {stage}")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Run one pmx NEQ manifest task.")
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--task-id", type=int, required=True)
    args = parser.parse_args(argv)

    try:
        run_task(args.manifest, args.task_id)
    except (GromacsError, FileNotFoundError, KeyError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
