from __future__ import annotations

"""Local MBAR analysis; no analysis dependencies are needed on Sherlock."""

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np


def analyze_phase(phase_dir: Path, temperature_k: float = 300.0) -> dict[str, float]:
    from pymbar import MBAR

    files = sorted(phase_dir.glob("state_*_energies.csv"))
    if not files:
        raise FileNotFoundError(f"No state energy CSV files in {phase_dir}")
    by_state: dict[int, list[list[float]]] = {}
    for path in files:
        with path.open() as handle:
            for row in csv.DictReader(handle):
                state = int(row["origin_state"])
                values = [float(row[key]) for key in row if key.startswith("u_")]
                by_state.setdefault(state, []).append(values)
    states = sorted(by_state)
    if states != list(range(len(states))):
        raise ValueError(f"Missing lambda windows; found states {states}")
    n_k = np.array([len(by_state[k]) for k in states], dtype=int)
    u_kn = np.concatenate(
        [np.asarray(by_state[k], dtype=float).T for k in states], axis=1
    )
    mbar = MBAR(u_kn, n_k)
    result = mbar.compute_free_energy_differences()
    delta = result["Delta_f"]
    uncertainty = result["dDelta_f"]
    rt_kj_mol = 0.00831446261815324 * temperature_k
    return {
        "delta_g_kj_mol": float(delta[0, -1] * rt_kj_mol),
        "uncertainty_kj_mol": float(uncertainty[0, -1] * rt_kj_mol),
        "samples": int(n_k.sum()),
        "minimum_samples_per_state": int(n_k.min()),
    }


def main() -> int:
    p = argparse.ArgumentParser(description="Analyze WT→V106A Jorgensen-style FEP")
    p.add_argument("--run-dir", type=Path, default=Path("results/analysis/fep_jorgensen/V106A"))
    args = p.parse_args()
    holo = analyze_phase(args.run_dir / "holo" / "windows")
    apo = analyze_phase(args.run_dir / "apo" / "windows")
    ddg = holo["delta_g_kj_mol"] - apo["delta_g_kj_mol"]
    unc = math.hypot(holo["uncertainty_kj_mol"], apo["uncertainty_kj_mol"])
    summary = {
        "reference": "WT",
        "mutation": "V106A",
        "wt_ddg_kj_mol": 0.0,
        "v106a_ddg_bind_kj_mol": ddg,
        "v106a_ddg_bind_kcal_mol": ddg / 4.184,
        "uncertainty_kj_mol": unc,
        "sign_convention": "positive means weaker DOR binding than WT",
        "holo": holo,
        "apo": apo,
    }
    destination = args.run_dir / "summary.json"
    destination.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
