from __future__ import annotations

import argparse
import csv
from pathlib import Path

from .worker import run_window


def read_task(manifest: Path, task_id: int) -> dict[str, str]:
    with manifest.open(newline="") as handle:
        for row in csv.DictReader(handle):
            if int(row["task_id"]) == task_id:
                return row
    raise IndexError(f"Task {task_id} not found in {manifest}")


def main() -> int:
    parser = argparse.ArgumentParser(description="Run one OpenMM task from a panel manifest")
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--task-id", type=int, required=True)
    args = parser.parse_args()
    row = read_task(args.manifest, args.task_id)
    run_window(
        phase_dir=Path(row["phase_dir"]),
        output_dir=Path(row["window_dir"]),
        state_index=int(row["state_index"]),
        temperature_k=float(row["temperature_k"]),
        timestep_fs=float(row["timestep_fs"]),
        collision_rate_per_ps=float(row["collision_rate_per_ps"]),
        equilibration_steps=int(row["equilibration_steps"]),
        production_steps=int(row["production_steps"]),
        energy_interval=int(row["energy_interval"]),
        checkpoint_interval=int(row["checkpoint_interval"]),
        platform_name=row["platform"],
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
