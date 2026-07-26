#!/usr/bin/env python3
"""Prepare NEQ directory layout, MDP files, and worker manifest for pmx FEP."""

from __future__ import annotations

import argparse
import csv
import json
import shutil
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.fep_pmx.config import (
    FEP_PMX_ROOT,
    NEQ_EQUIL_NS,
    NEQ_EQUIL_SNAPSHOT_START_PS,
    NEQ_SNAPSHOTS_DEFAULT,
    P0_LEGS,
    switch_bundle_ranges,
    switch_ps_for_leg,
)
from scripts.fep_pmx.mdp_utils import render_em_mdp, render_npt_eq_mdp, render_nonequil_mdp


def _leg_dir(leg_id: str, phase: str, replicate: int) -> Path:
    return FEP_PMX_ROOT / "legs" / leg_id / phase / f"rep_{replicate:02d}"


def _build_dir(leg_id: str, phase: str, replicate: int) -> Path:
    return _leg_dir(leg_id, phase, replicate) / "gromacs_build"


def _neq_dir(leg_id: str, phase: str, replicate: int) -> Path:
    return _leg_dir(leg_id, phase, replicate) / "neq"


def _snapshot_times_ps(n_snapshots: int) -> list[float]:
    """Evenly spaced snapshot times during equil production window."""
    start_ps = NEQ_EQUIL_SNAPSHOT_START_PS
    end_ps = NEQ_EQUIL_NS * 1000.0
    if n_snapshots < 1:
        raise ValueError("n_snapshots must be >= 1")
    if end_ps <= start_ps:
        raise ValueError("equil window too short for snapshot extraction")
    if n_snapshots == 1:
        return [start_ps]
    step = (end_ps - start_ps) / float(n_snapshots - 1)
    return [start_ps + i * step for i in range(n_snapshots)]


def _copy_topology_bundle(build_dir: Path, neq_dir: Path, top_path: Path) -> None:
    """Copy local topology dependencies from gromacs_build into neq/.

    pmx gentop writes chain itps (pmx_topol_*.itp) referenced transitively from
    system.top; copying every *.itp from the build dir avoids missing nested includes.
    """
    copied: set[str] = set()
    for src in sorted(build_dir.glob("*.itp")):
        shutil.copy2(src, neq_dir / src.name)
        copied.add(src.name)

    for line in top_path.read_text().splitlines():
        stripped = line.strip()
        if not stripped.startswith("#include"):
            continue
        parts = stripped.split()
        if len(parts) < 2:
            continue
        include_name = parts[1].strip('"')
        if "/" in include_name or include_name.startswith("."):
            continue
        if include_name in copied:
            continue
        src = build_dir / include_name
        if not src.is_file():
            raise FileNotFoundError(
                f"system.top references missing include {include_name} under {build_dir}"
            )
        shutil.copy2(src, neq_dir / include_name)
        copied.add(include_name)


def prepare_neq(
    leg_id: str,
    *,
    phase: str,
    replicate: int = 1,
    n_snapshots: int = NEQ_SNAPSHOTS_DEFAULT,
    force: bool = False,
) -> Path:
    build = _build_dir(leg_id, phase, replicate)
    gro = build / "system.gro"
    top = build / "system.top"
    if not gro.is_file() or not top.is_file():
        raise FileNotFoundError(f"Missing solvated system under {build}")

    neq = _neq_dir(leg_id, phase, replicate)
    if neq.is_dir() and (neq / "neq_manifest.csv").is_file() and not force:
        neq.mkdir(parents=True, exist_ok=True)
        shutil.copy2(gro, neq / "system.gro")
        shutil.copy2(top, neq / "system.top")
        _copy_topology_bundle(build, neq, top)
        return neq / "neq_manifest.csv"

    # force: refresh mdps/manifests only — never delete em/equil/extract/switch outputs.
    neq.mkdir(parents=True, exist_ok=True)

    mdp_dir = neq / "mdp"
    mdp_dir.mkdir(exist_ok=True)
    render_em_mdp(mdp_dir / "em.mdp")
    render_npt_eq_mdp(output=mdp_dir / "npt_eq_lambda0.mdp", init_lambda=0.0)
    render_npt_eq_mdp(output=mdp_dir / "npt_eq_lambda1.mdp", init_lambda=1.0)
    switch_ps = switch_ps_for_leg(leg_id)
    render_nonequil_mdp(output=mdp_dir / "nonequil_fwd.mdp", init_lambda=0.0, switch_ps=switch_ps)
    render_nonequil_mdp(output=mdp_dir / "nonequil_rev.mdp", init_lambda=1.0, switch_ps=switch_ps)

    for sub in ("em", "eq_lambda0", "eq_lambda1", "snapshots/lambda0", "snapshots/lambda1", "switches"):
        (neq / sub).mkdir(parents=True, exist_ok=True)

    shutil.copy2(gro, neq / "system.gro")
    shutil.copy2(top, neq / "system.top")
    _copy_topology_bundle(build, neq, top)

    snapshot_times = _snapshot_times_ps(n_snapshots)
    manifest_rows: list[dict[str, str | int | float]] = []
    task_id = 0

    manifest_rows.append(
        {
            "task_id": task_id,
            "leg_id": leg_id,
            "phase": phase,
            "replicate": replicate,
            "stage": "em",
            "lambda_state": "",
            "direction": "",
            "snapshot_index": "",
            "snapshot_index_end": "",
            "snapshot_time_ps": "",
            "switch_ps": "",
            "run_dir": "em",
        }
    )
    task_id += 1

    for lambda_state in (0, 1):
        manifest_rows.append(
            {
                "task_id": task_id,
                "leg_id": leg_id,
                "phase": phase,
                "replicate": replicate,
                "stage": "equil",
                "lambda_state": lambda_state,
                "direction": "",
                "snapshot_index": "",
                "snapshot_index_end": "",
                "snapshot_time_ps": "",
                "switch_ps": "",
                "run_dir": f"eq_lambda{lambda_state}",
            }
        )
        task_id += 1

        manifest_rows.append(
            {
                "task_id": task_id,
                "leg_id": leg_id,
                "phase": phase,
                "replicate": replicate,
                "stage": "extract",
                "lambda_state": lambda_state,
                "direction": "",
                "snapshot_index": "",
                "snapshot_index_end": "",
                "snapshot_time_ps": "",
                "switch_ps": "",
                "run_dir": f"snapshots/lambda{lambda_state}",
            }
        )
        task_id += 1

    direction_by_lambda = {0: "fwd", 1: "rev"}
    for lambda_state in (0, 1):
        direction = direction_by_lambda[lambda_state]
        for start_idx, end_idx in switch_bundle_ranges(n_snapshots, leg_id):
            manifest_rows.append(
                {
                    "task_id": task_id,
                    "leg_id": leg_id,
                    "phase": phase,
                    "replicate": replicate,
                    "stage": "switch",
                    "lambda_state": lambda_state,
                    "direction": direction,
                    "snapshot_index": start_idx,
                    "snapshot_index_end": end_idx,
                    "snapshot_time_ps": f"{snapshot_times[start_idx]:.3f}",
                    "switch_ps": switch_ps,
                    "run_dir": f"switches/{direction}_bundle{start_idx:03d}",
                }
            )
            task_id += 1
            for idx in range(start_idx, end_idx + 1):
                (neq / "switches" / f"{direction}_{idx:03d}").mkdir(parents=True, exist_ok=True)

    manifest_path = neq / "neq_manifest.csv"
    fieldnames = list(manifest_rows[0].keys())
    with manifest_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(manifest_rows)

    meta = {
        "leg_id": leg_id,
        "phase": phase,
        "replicate": replicate,
        "n_snapshots": n_snapshots,
        "switch_ps": switch_ps,
        "equil_ns": NEQ_EQUIL_NS,
        "snapshot_times_ps": snapshot_times,
        "system_gro": str(gro),
        "system_top": str(top),
        "manifest": str(manifest_path),
        "n_tasks": len(manifest_rows),
    }
    (neq / "neq_prepare.json").write_text(json.dumps(meta, indent=2) + "\n")
    return manifest_path


def refresh_switch_schedule(
    leg_id: str,
    *,
    phase: str,
    replicate: int,
    n_snapshots: int,
) -> Path:
    """Update switch task schedule without rebuilding em/equil/extract (safe mid-pipeline)."""
    neq = _neq_dir(leg_id, phase, replicate)
    manifest_path = neq / "neq_manifest.csv"
    if not manifest_path.is_file():
        raise FileNotFoundError(f"Missing NEQ manifest: {manifest_path}")

    with manifest_path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        fieldnames = reader.fieldnames
        if not fieldnames:
            raise ValueError(f"Empty manifest: {manifest_path}")
        rows = list(reader)

    prefix = [row for row in rows if row["stage"] != "switch"]
    if not prefix:
        raise ValueError(f"No em/equil/extract rows in {manifest_path}")

    switch_ps = switch_ps_for_leg(leg_id)
    snapshot_times = _snapshot_times_ps(n_snapshots)
    task_id = max(int(row["task_id"]) for row in prefix) + 1
    direction_by_lambda = {0: "fwd", 1: "rev"}

    for lambda_state in (0, 1):
        direction = direction_by_lambda[lambda_state]
        snap_dir = neq / f"snapshots/lambda{lambda_state}"
        snap_dir.mkdir(parents=True, exist_ok=True)
        for start_idx, end_idx in switch_bundle_ranges(n_snapshots, leg_id):
            for idx in range(start_idx, end_idx + 1):
                (neq / "switches" / f"{direction}_{idx:03d}").mkdir(parents=True, exist_ok=True)
            prefix.append(
                {
                    "task_id": task_id,
                    "leg_id": leg_id,
                    "phase": phase,
                    "replicate": replicate,
                    "stage": "switch",
                    "lambda_state": lambda_state,
                    "direction": direction,
                    "snapshot_index": start_idx,
                    "snapshot_index_end": end_idx,
                    "snapshot_time_ps": f"{snapshot_times[start_idx]:.3f}",
                    "switch_ps": switch_ps,
                    "run_dir": f"switches/{direction}_bundle{start_idx:03d}",
                }
            )
            task_id += 1

    with manifest_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(prefix)

    meta_path = neq / "neq_prepare.json"
    meta = json.loads(meta_path.read_text()) if meta_path.is_file() else {}
    meta.update(
        {
            "leg_id": leg_id,
            "phase": phase,
            "replicate": replicate,
            "n_snapshots": n_snapshots,
            "switch_ps": switch_ps,
            "equil_ns": NEQ_EQUIL_NS,
            "snapshot_times_ps": snapshot_times,
            "manifest": str(manifest_path),
            "n_tasks": len(prefix),
        }
    )
    meta_path.write_text(json.dumps(meta, indent=2) + "\n")
    return manifest_path


def rebuild_panel_manifest(
    *,
    legs: tuple[str, ...],
    phases: tuple[str, ...],
    replicates: range,
    output: Path,
) -> Path:
    """Rebuild panel manifest from existing per-leg neq_manifest.csv files."""
    rows: list[dict[str, str]] = []
    offset = 0
    for leg_id in legs:
        for phase in phases:
            for replicate in replicates:
                manifest = _neq_dir(leg_id, phase, replicate) / "neq_manifest.csv"
                if not manifest.is_file():
                    raise FileNotFoundError(f"Missing per-leg manifest: {manifest}")
                with manifest.open(newline="") as handle:
                    reader = csv.DictReader(handle)
                    for row in reader:
                        row = dict(row)
                        row["panel_task_id"] = str(offset + int(row["task_id"]))
                        rows.append(row)
                offset = len(rows)

    if not rows:
        raise ValueError("No manifest rows found")

    fieldnames = ["panel_task_id", *[k for k in rows[0].keys() if k != "panel_task_id"]]
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    return output


def build_panel_manifest(
    *,
    legs: tuple[str, ...],
    phases: tuple[str, ...],
    replicates: range,
    n_snapshots: int,
    output: Path,
    force: bool = False,
) -> Path:
    per_leg_manifests: list[Path] = []
    for leg_id in legs:
        for phase in phases:
            for replicate in replicates:
                per_leg_manifests.append(
                    prepare_neq(
                        leg_id,
                        phase=phase,
                        replicate=replicate,
                        n_snapshots=n_snapshots,
                        force=force,
                    )
                )

    rows: list[dict[str, str]] = []
    offset = 0
    for manifest in per_leg_manifests:
        with manifest.open(newline="") as handle:
            reader = csv.DictReader(handle)
            for row in reader:
                row = dict(row)
                row["panel_task_id"] = str(offset + int(row["task_id"]))
                rows.append(row)
        offset = len(rows)

    fieldnames = ["panel_task_id", *[k for k in rows[0].keys() if k != "panel_task_id"]]
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    return output


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Prepare pmx NEQ manifests and MDP inputs.")
    parser.add_argument("--leg", default=None, help="Single leg id (omit for P0 panel manifest)")
    parser.add_argument("--phase", choices=("holo", "apo"), default="holo")
    parser.add_argument("--replicate", type=int, default=1)
    parser.add_argument("--n-snapshots", type=int, default=NEQ_SNAPSHOTS_DEFAULT)
    parser.add_argument("--force", action="store_true")
    parser.add_argument(
        "--refresh-switch-only",
        action="store_true",
        help="Update switch schedule + snapshot times only (keep em/equil/extract; safe while jobs run)",
    )
    parser.add_argument(
        "--rebuild-panel-only",
        action="store_true",
        help="Rebuild panel CSV from existing per-leg manifests (after refresh-switch-only)",
    )
    parser.add_argument(
        "--panel-manifest",
        type=Path,
        default=FEP_PMX_ROOT / "neq_panel_manifest.csv",
        help="Write combined manifest for P0 legs when --leg is omitted",
    )
    parser.add_argument("--replicates", type=int, default=1, help="Panel mode: reps 1..N")
    args = parser.parse_args(argv)

    if args.rebuild_panel_only:
        manifest = rebuild_panel_manifest(
            legs=P0_LEGS,
            phases=("holo", "apo"),
            replicates=range(1, args.replicates + 1),
            output=args.panel_manifest,
        )
        print(f"Rebuilt panel NEQ manifest: {manifest}")
        return 0

    if args.refresh_switch_only:
        if args.leg:
            refresh_switch_schedule(
                args.leg,
                phase=args.phase,
                replicate=args.replicate,
                n_snapshots=args.n_snapshots,
            )
            print(f"Refreshed switch schedule: {args.leg} {args.phase} rep{args.replicate}")
            return 0
        for leg_id in P0_LEGS:
            for phase in ("holo", "apo"):
                for replicate in range(1, args.replicates + 1):
                    refresh_switch_schedule(
                        leg_id,
                        phase=phase,
                        replicate=replicate,
                        n_snapshots=args.n_snapshots,
                    )
        print(
            f"Refreshed switch schedule for P0 panel "
            f"({args.n_snapshots} snapshots / {args.replicates} reps)"
        )
        return 0

    if args.leg:
        manifest = prepare_neq(
            args.leg,
            phase=args.phase,
            replicate=args.replicate,
            n_snapshots=args.n_snapshots,
            force=args.force,
        )
        print(f"Wrote NEQ manifest: {manifest}")
        return 0

    manifest = build_panel_manifest(
        legs=P0_LEGS,
        phases=("holo", "apo"),
        replicates=range(1, args.replicates + 1),
        n_snapshots=args.n_snapshots,
        output=args.panel_manifest,
        force=args.force,
    )
    print(f"Wrote panel NEQ manifest: {manifest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
