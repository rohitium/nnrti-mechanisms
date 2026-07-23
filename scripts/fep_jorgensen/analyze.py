"""Local MBAR analysis; no analysis dependencies are needed on Sherlock."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np

from .mutations import MANUSCRIPT_PLANS, canonical_label, safe_label


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


def analyze_leg(run_dir: Path, temperature_k: float | None = None) -> dict:
    config_path = run_dir / "config.json"
    config = json.loads(config_path.read_text()) if config_path.exists() else {}
    temperature = float(temperature_k or config.get("temperature_k", 300.0))
    holo = analyze_phase(run_dir / "holo" / "windows", temperature)
    apo = analyze_phase(run_dir / "apo" / "windows", temperature)
    ddg = holo["delta_g_kj_mol"] - apo["delta_g_kj_mol"]
    uncertainty = math.hypot(
        holo["uncertainty_kj_mol"], apo["uncertainty_kj_mol"]
    )
    summary = {
        "leg_id": config.get("leg_id", run_dir.name),
        "start_label": config.get("start_label"),
        "end_label": config.get("end_label"),
        "mutation": config.get("mutation"),
        "ddg_bind_kj_mol": ddg,
        "ddg_bind_kcal_mol": ddg / 4.184,
        "uncertainty_kj_mol": uncertainty,
        "uncertainty_kcal_mol": uncertainty / 4.184,
        "sign_convention": "positive means weaker DOR binding at the end state",
        "holo": holo,
        "apo": apo,
    }
    (run_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    return summary


def analyze_target(target: str, output_dir: Path) -> dict:
    target = canonical_label(target)
    try:
        plan = MANUSCRIPT_PLANS[target]
    except KeyError as exc:
        raise ValueError(f"Target {target} is not in the manuscript panel") from exc
    leg_summaries = []
    for leg in plan.legs:
        run_dir = output_dir / "legs" / leg.leg_id
        summary_path = run_dir / "summary.json"
        summary = (
            json.loads(summary_path.read_text())
            if summary_path.exists()
            else analyze_leg(run_dir)
        )
        leg_summaries.append(summary)
    ddg = sum(float(leg["ddg_bind_kj_mol"]) for leg in leg_summaries)
    uncertainty = math.sqrt(
        sum(float(leg["uncertainty_kj_mol"]) ** 2 for leg in leg_summaries)
    )
    result = {
        "reference": "WT",
        "target": target,
        "ddg_bind_kj_mol": ddg,
        "ddg_bind_kcal_mol": ddg / 4.184,
        "uncertainty_kj_mol": uncertainty,
        "uncertainty_kcal_mol": uncertainty / 4.184,
        "sign_convention": "positive means weaker DOR binding than WT",
        "legs": leg_summaries,
    }
    destination = output_dir / "targets" / safe_label(target) / "summary.json"
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(result, indent=2) + "\n")
    return result


def main() -> int:
    p = argparse.ArgumentParser(description="Analyze mutation-agnostic Jorgensen-style FEP")
    group = p.add_mutually_exclusive_group(required=True)
    group.add_argument("--leg-dir", "--run-dir", dest="leg_dir", type=Path)
    group.add_argument("--target")
    group.add_argument("--all-targets", action="store_true")
    p.add_argument("--output-dir", type=Path, default=Path("results/analysis/fep_jorgensen"))
    args = p.parse_args()
    if args.leg_dir:
        result = analyze_leg(args.leg_dir)
    elif args.target:
        result = analyze_target(args.target, args.output_dir)
    else:
        result = {
            target: analyze_target(target, args.output_dir)
            for target in MANUSCRIPT_PLANS
        }
        table = args.output_dir / "manuscript_panel_summary.csv"
        with table.open("w", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=["target", "ddg_bind_kj_mol", "ddg_bind_kcal_mol",
                            "uncertainty_kj_mol", "uncertainty_kcal_mol"],
            )
            writer.writeheader()
            for summary in result.values():
                writer.writerow({key: summary[key] for key in writer.fieldnames})
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
