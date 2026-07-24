"""Convergence diagnostics for fixed-lambda FEP window outputs."""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np


def _window_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def diagnose_windows(
    windows_dir: Path,
    *,
    temperature_k: float = 300.0,
    target_samples: int = 1000,
    drift_fraction: float = 0.05,
) -> dict:
    """Summarize sample counts and energy drift for one phase."""
    files = sorted(windows_dir.glob("state_*_energies.csv"))
    if not files:
        raise FileNotFoundError(f"No state energy CSV files in {windows_dir}")

    rt_kj_mol = 0.00831446261815324 * temperature_k
    windows: list[dict] = []
    for path in files:
        state = int(path.stem.split("_")[1])
        rows = _window_rows(path)
        sample_count = len(rows)
        status = "ok"
        if sample_count < target_samples:
            status = "partial" if sample_count else "missing"
        drift_kj_mol = None
        if sample_count >= 2:
            key = f"u_{state}"
            own = np.asarray([float(row[key]) for row in rows], dtype=float)
            block = max(1, int(sample_count * drift_fraction))
            drift_kj_mol = float((own[-block:].mean() - own[:block].mean()) * rt_kj_mol)
            if drift_kj_mol is not None and abs(drift_kj_mol) > 5.0 and status == "ok":
                status = "drifting"
        windows.append(
            {
                "state_index": state,
                "samples": sample_count,
                "target_samples": target_samples,
                "drift_kj_mol": drift_kj_mol,
                "status": status,
            }
        )

    complete = sum(1 for row in windows if row["status"] == "ok")
    return {
        "windows_dir": str(windows_dir),
        "n_windows": len(windows),
        "complete_windows": complete,
        "minimum_samples": min(row["samples"] for row in windows),
        "maximum_drift_kj_mol": max(
            (abs(row["drift_kj_mol"]) for row in windows if row["drift_kj_mol"] is not None),
            default=0.0,
        ),
        "windows": windows,
        "ready_for_mbar": complete == len(windows),
    }


def diagnose_phase(phase_dir: Path, **kwargs) -> dict:
    return diagnose_windows(phase_dir / "windows", **kwargs)


def diagnose_leg(run_dir: Path, **kwargs) -> dict:
    phases: dict[str, dict] = {}
    for phase in ("holo", "apo"):
        windows = run_dir / phase / "windows"
        if windows.is_dir() and any(windows.glob("state_*_energies.csv")):
            phases[phase] = diagnose_phase(run_dir / phase, **kwargs)
    return {"leg_id": run_dir.name, "phases": phases}
