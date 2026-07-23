"""Local MBAR analysis for holo-only mutation FEP legs."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np

from .mutations import MANUSCRIPT_PLANS, canonical_label, safe_label


def analyze_phase(phase_dir: Path, temperature_k: float = 300.0) -> dict[str, float]:
    windows = phase_dir / "windows"
    if windows.is_dir() and any(windows.glob("state_*_energies.csv")):
        return _analyze_window_energies(windows, temperature_k)
    reporter_path = phase_dir / "multistate.nc"
    if reporter_path.is_file():
        return _analyze_multistate_reporter(reporter_path, temperature_k)
    raise FileNotFoundError(
        f"No MBAR inputs in {phase_dir}; expected windows/state_*_energies.csv or multistate.nc"
    )


def _analyze_window_energies(windows_dir: Path, temperature_k: float) -> dict[str, float]:
    from pymbar import MBAR

    files = sorted(windows_dir.glob("state_*_energies.csv"))
    if not files:
        raise FileNotFoundError(f"No state energy CSV files in {windows_dir}")
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


def _analyze_multistate_reporter(reporter_path: Path, temperature_k: float) -> dict[str, float]:
    from openmm import unit
    from openmmtools.constants import kB
    from openmmtools.multistate import MultiStateReporter, MultiStateSamplerAnalyzer

    reporter = MultiStateReporter(str(reporter_path), open_mode="r")
    analyzer = MultiStateSamplerAnalyzer(reporter)
    delta_f, d_delta_f = analyzer.get_free_energy()
    kT = (kB * temperature_k * unit.kelvin).value_in_unit(unit.kilojoule_per_mole)
    reporter.close()
    return {
        "delta_g_kj_mol": float(delta_f[0, -1] * kT),
        "uncertainty_kj_mol": float(d_delta_f[0, -1] * kT),
        "samples": -1,
        "minimum_samples_per_state": -1,
    }


def analyze_leg(run_dir: Path, temperature_k: float | None = None) -> dict:
    config_path = run_dir / "config.json"
    config = json.loads(config_path.read_text()) if config_path.exists() else {}
    temperature = float(temperature_k or config.get("temperature_k", 300.0))
    holo = analyze_phase(run_dir / "holo", temperature)
    holo_schedule = json.loads((run_dir / "holo" / "schedule.json").read_text())
    backend = holo_schedule.get("prepare_backend", "scaling")
    holo_strategy = holo_schedule.get("alchemical_plan", {}).get(
        "strategy", "annihilate_wt_sidechain"
    )
    delta_kj = holo["delta_g_kj_mol"]
    uncertainty_kj = holo["uncertainty_kj_mol"]
    if backend == "scaling" and holo_strategy == "annihilate_mutant_sidechain":
        delta_kj *= -1.0
    summary = {
        "leg_id": config.get("leg_id", run_dir.name),
        "start_label": config.get("start_label"),
        "end_label": config.get("end_label"),
        "mutation": config.get("mutation"),
        "delta_g_mutation_kj_mol": delta_kj,
        "delta_g_mutation_kcal_mol": delta_kj / 4.184,
        "uncertainty_kj_mol": uncertainty_kj,
        "uncertainty_kcal_mol": uncertainty_kj / 4.184,
        "thermodynamic_cycle": "protein-side-chain mutation in inhibitor-bound complex",
        "sign_convention": "positive means the end-state mutation is higher in free energy",
        "prepare_backend": backend,
        "sampling_strategy": holo_strategy if backend == "scaling" else "perses-default",
        "holo": holo,
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
    delta_g = sum(float(leg["delta_g_mutation_kj_mol"]) for leg in leg_summaries)
    uncertainty = math.sqrt(
        sum(float(leg["uncertainty_kj_mol"]) ** 2 for leg in leg_summaries)
    )
    result = {
        "reference": "WT",
        "target": target,
        "delta_g_mutation_kj_mol": delta_g,
        "delta_g_mutation_kcal_mol": delta_g / 4.184,
        "uncertainty_kj_mol": uncertainty,
        "uncertainty_kcal_mol": uncertainty / 4.184,
        "thermodynamic_cycle": "protein-side-chain mutation in inhibitor-bound complex",
        "sign_convention": "positive means the target mutation costs more free energy than WT",
        "legs": leg_summaries,
    }
    destination = output_dir / "targets" / safe_label(target) / "summary.json"
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(result, indent=2) + "\n")
    return result


def normalize_to_reference(
    summaries: dict[str, dict], reference: str = "WT"
) -> list[dict[str, float | str]]:
    reference = canonical_label(reference)
    if reference in summaries:
        reference_dg = float(summaries[reference]["delta_g_mutation_kcal_mol"])
        reference_sigma = float(summaries[reference]["uncertainty_kcal_mol"])
    elif reference == "WT":
        reference_dg = 0.0
        reference_sigma = 0.0
    else:
        raise ValueError(f"Reference target {reference} is missing from summaries")
    rows = []
    for target, summary in sorted(summaries.items()):
        delta_g = float(summary["delta_g_mutation_kcal_mol"])
        sigma = float(summary["uncertainty_kcal_mol"])
        is_reference = canonical_label(target) == reference
        rows.append(
            {
                "reference": reference,
                "target": target,
                "delta_delta_g_kcal_mol": 0.0 if is_reference else delta_g - reference_dg,
                "uncertainty_kcal_mol": 0.0 if is_reference else math.hypot(sigma, reference_sigma),
            }
        )
    return rows


def main() -> int:
    p = argparse.ArgumentParser(
        description="Analyze holo-only Jorgensen-inspired mutation FEP"
    )
    group = p.add_mutually_exclusive_group(required=True)
    group.add_argument("--leg-dir", "--run-dir", dest="leg_dir", type=Path)
    group.add_argument("--target")
    group.add_argument("--all-targets", action="store_true")
    p.add_argument("--output-dir", type=Path, default=Path("results/analysis/fep_jorgensen"))
    p.add_argument("--reference", default="WT", help="Reference system for relative normalization")
    p.add_argument("--relative-table", type=Path, help="Optional CSV for reference-normalized values")
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
                fieldnames=[
                    "target", "delta_g_mutation_kj_mol", "delta_g_mutation_kcal_mol",
                    "uncertainty_kj_mol", "uncertainty_kcal_mol",
                ],
            )
            writer.writeheader()
            for summary in result.values():
                writer.writerow({key: summary[key] for key in writer.fieldnames})
        relative_rows = normalize_to_reference(result, reference=args.reference)
        relative_table = args.relative_table or args.output_dir / "manuscript_panel_relative.csv"
        with relative_table.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(relative_rows[0]))
            writer.writeheader()
            writer.writerows(relative_rows)
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
