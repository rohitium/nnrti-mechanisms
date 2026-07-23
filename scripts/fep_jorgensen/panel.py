from __future__ import annotations

import argparse
import csv
from pathlib import Path
import shlex

from .config import FEPConfig
from .mutations import MANUSCRIPT_TARGETS, unique_manuscript_legs


def preparation_commands(output_dir: Path, replicate: int = 1) -> list[str]:
    commands = []
    for leg in unique_manuscript_legs():
        commands.append(
            shlex.join(
                [
                    "python", "-m", "scripts.fep_jorgensen.prepare",
                    "--mutation", leg.mutation,
                    "--start-label", leg.start_label,
                    "--end-label", leg.end_label,
                    "--input-complex-pdb", str(leg.input_complex_pdb(replicate)),
                    "--output-dir", str(output_dir),
                ]
            )
        )
    return commands


def write_worker_manifest(
    destination: Path,
    output_dir: Path,
    config: FEPConfig | None = None,
) -> int:
    settings = config or FEPConfig(output_dir=output_dir)
    fields = [
        "task_id", "leg_id", "start_label", "end_label", "mutation", "phase",
        "state_index", "phase_dir", "window_dir", "temperature_k", "timestep_fs",
        "collision_rate_per_ps", "equilibration_steps", "production_steps",
        "energy_interval", "checkpoint_interval", "platform",
    ]
    destination.parent.mkdir(parents=True, exist_ok=True)
    task_id = 0
    with destination.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for leg in unique_manuscript_legs():
            leg_dir = output_dir / "legs" / leg.leg_id
            for state_index in range(len(settings.lambda_schedule.values)):
                writer.writerow(
                    {
                        "task_id": task_id,
                        "leg_id": leg.leg_id,
                        "start_label": leg.start_label,
                        "end_label": leg.end_label,
                        "mutation": leg.mutation,
                        "phase": "holo",
                        "state_index": state_index,
                        "phase_dir": leg_dir / "holo",
                        "window_dir": leg_dir / "holo" / "windows",
                        "temperature_k": settings.temperature_k,
                        "timestep_fs": settings.timestep_fs,
                        "collision_rate_per_ps": settings.collision_rate_per_ps,
                        "equilibration_steps": settings.equilibration_steps,
                        "production_steps": settings.production_steps,
                        "energy_interval": settings.energy_interval,
                        "checkpoint_interval": settings.checkpoint_interval,
                        "platform": settings.platform,
                    }
                )
                task_id += 1
    return task_id


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Plan the merged Jorgensen-inspired manuscript FEP panel"
    )
    parser.add_argument("--output-dir", type=Path, default=FEPConfig().output_dir)
    parser.add_argument("--replicate", type=int, default=1)
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--preparation-script", type=Path)
    args = parser.parse_args()
    manifest = args.manifest or args.output_dir / "worker_manifest.csv"
    count = write_worker_manifest(manifest, args.output_dir)
    preparation_script = args.preparation_script or args.output_dir / "prepare_all.sh"
    preparation_script.parent.mkdir(parents=True, exist_ok=True)
    preparation_script.write_text(
        "#!/bin/bash\nset -euo pipefail\n\n"
        + "\n".join(preparation_commands(args.output_dir, args.replicate))
        + "\n"
    )
    preparation_script.chmod(0o755)
    print(f"Targets: {len(MANUSCRIPT_TARGETS)}")
    print(f"Unique alchemical legs: {len(unique_manuscript_legs())}")
    print(f"OpenMM holo worker tasks: {count}")
    print(f"Preparation script: {preparation_script}")
    print(f"Worker manifest: {manifest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
