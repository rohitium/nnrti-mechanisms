from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
import pytest

from scripts.fep_jorgensen.analyze import analyze_phase, analyze_target
from scripts.fep_jorgensen.config import FEPConfig, LambdaSchedule
from scripts.fep_jorgensen.mutations import (
    MANUSCRIPT_PLANS,
    MANUSCRIPT_TARGETS,
    Mutation,
    MutationLeg,
    unique_manuscript_legs,
)
from scripts.fep_jorgensen.panel import write_worker_manifest
from scripts.fep_jorgensen.worker import _perses_default_parameters, run_window


def _write_toy_phase(path: Path) -> None:
    from openmm import CustomExternalForce, System, XmlSerializer

    path.mkdir(parents=True)
    system = System()
    system.addParticle(39.9)
    force = CustomExternalForce("lambda_sterics_core")
    force.addGlobalParameter("lambda_sterics_core", 0.0)
    force.addParticle(0, [])
    system.addForce(force)
    (path / "hybrid_system.xml").write_text(XmlSerializer.serialize(system))
    (path / "hybrid_topology.pdb").write_text(
        "ATOM      1  AR  ARG A   1       0.000   0.000   0.000  1.00  0.00          Ar\n"
        "TER\nEND\n"
    )
    (path / "schedule.json").write_text(json.dumps({"lambda_values": [0.0, 0.5, 1.0]}))


def _worker_arguments(phase: Path, windows: Path) -> dict:
    return {
        "phase_dir": phase,
        "output_dir": windows,
        "state_index": 0,
        "temperature_k": 300.0,
        "timestep_fs": 1.0,
        "collision_rate_per_ps": 1.0,
        "equilibration_steps": 1,
        "production_steps": 4,
        "energy_interval": 1,
        "checkpoint_interval": 2,
        "platform_name": "CPU",
    }


def test_configuration_validation() -> None:
    LambdaSchedule((0.0, 0.5, 1.0)).validate()
    with pytest.raises(ValueError):
        LambdaSchedule((0.0, 0.5, 0.5, 1.0)).validate()
    with pytest.raises(ValueError):
        FEPConfig(mutation="Y181C").validate()
    config = FEPConfig.for_leg(MANUSCRIPT_PLANS["Y181C"].legs[0])
    config.validate()
    assert (config.residue_id, config.wt_residue, config.mutant_residue) == (
        "181", "TYR", "CYS"
    )
    assert config.run_dir.name == "wt_to_Y181C"


def test_complete_manuscript_panel_has_continuous_single_residue_legs() -> None:
    assert len(MANUSCRIPT_TARGETS) == 19
    assert len(unique_manuscript_legs()) == 19
    for target, plan in MANUSCRIPT_PLANS.items():
        assert plan.legs[0].start_label == "WT"
        assert plan.legs[-1].end_label == target
        for leg in plan.legs:
            Mutation.parse(leg.mutation)
        for first, second in zip(plan.legs, plan.legs[1:]):
            assert first.end_label == second.start_label
    with pytest.raises(ValueError, match="does not add exactly"):
        MutationLeg("WT", "Y181C", "V106A")


def test_panel_manifest_has_every_phase_and_lambda_state(tmp_path: Path) -> None:
    manifest = tmp_path / "worker_manifest.csv"
    config = FEPConfig(
        output_dir=tmp_path,
        lambda_schedule=LambdaSchedule((0.0, 0.5, 1.0)),
    )
    count = write_worker_manifest(manifest, tmp_path, config)
    rows = list(csv.DictReader(manifest.open()))
    assert count == 19 * 2 * 3
    assert len(rows) == count
    assert {row["phase"] for row in rows} == {"holo", "apo"}
    assert {int(row["state_index"]) for row in rows} == {0, 1, 2}


def test_double_mutant_analysis_sums_sequential_legs(tmp_path: Path) -> None:
    plan = MANUSCRIPT_PLANS["V106A+L234I"]
    values = ((2.0, 0.3), (4.0, 0.4))
    for leg, (ddg, uncertainty) in zip(plan.legs, values):
        run_dir = tmp_path / "legs" / leg.leg_id
        run_dir.mkdir(parents=True)
        (run_dir / "summary.json").write_text(
            json.dumps(
                {
                    "leg_id": leg.leg_id,
                    "start_label": leg.start_label,
                    "end_label": leg.end_label,
                    "mutation": leg.mutation,
                    "ddg_bind_kj_mol": ddg,
                    "uncertainty_kj_mol": uncertainty,
                }
            )
        )
    result = analyze_target("V106A_L234I", tmp_path)
    assert result["ddg_bind_kj_mol"] == pytest.approx(6.0)
    assert result["uncertainty_kj_mol"] == pytest.approx(0.5)


@pytest.mark.parametrize(
    ("master", "expected"),
    [(0.0, (0.0, 0.0)), (0.25, (0.5, 0.0)), (0.5, (1.0, 0.0)),
     (0.75, (1.0, 0.5)), (1.0, (1.0, 1.0))],
)
def test_perses_default_staging(master: float, expected: tuple[float, float]) -> None:
    values = _perses_default_parameters(master)
    assert values["lambda_sterics_insert"] == expected[0]
    assert values["lambda_electrostatics_delete"] == expected[0]
    assert values["lambda_sterics_delete"] == expected[1]
    assert values["lambda_electrostatics_insert"] == expected[1]


def test_openmm_worker_and_mbar_round_trip(tmp_path: Path) -> None:
    phase, windows = tmp_path / "phase", tmp_path / "windows"
    _write_toy_phase(phase)
    for state in range(3):
        arguments = _worker_arguments(phase, windows)
        arguments["state_index"] = state
        run_window(**arguments)
    for state in range(3):
        rows = list(csv.DictReader((windows / f"state_{state:02d}_energies.csv").open()))
        assert len(rows) == 4
        assert [float(rows[0][f"u_{i}"]) for i in range(3)] == pytest.approx(
            [0.0, 0.5 / (0.00831446261815324 * 300), 1.0 / (0.00831446261815324 * 300)]
        )
    result = analyze_phase(windows, temperature_k=300.0)
    assert result["samples"] == 12
    assert result["minimum_samples_per_state"] == 4
    assert result["delta_g_kj_mol"] == pytest.approx(1.0, abs=1e-6)


def test_completed_window_is_not_duplicated_on_resume(tmp_path: Path) -> None:
    phase, windows = tmp_path / "phase", tmp_path / "windows"
    _write_toy_phase(phase)
    arguments = _worker_arguments(phase, windows)
    run_window(**arguments)
    run_window(**arguments)
    rows = list(csv.DictReader((windows / "state_00_energies.csv").open()))
    assert len(rows) == 4
    assert [int(row["sample"]) for row in rows] == [0, 1, 2, 3]


def test_checkpoint_interval_must_align_with_energy_interval(tmp_path: Path) -> None:
    phase = tmp_path / "phase"
    _write_toy_phase(phase)
    arguments = _worker_arguments(phase, tmp_path / "windows")
    arguments.update(energy_interval=2, checkpoint_interval=3)
    with pytest.raises(ValueError, match="checkpoint_interval"):
        run_window(**arguments)


def test_production_steps_must_align_with_energy_interval(tmp_path: Path) -> None:
    phase = tmp_path / "phase"
    _write_toy_phase(phase)
    arguments = _worker_arguments(phase, tmp_path / "windows")
    arguments.update(production_steps=3, energy_interval=2, checkpoint_interval=2)
    with pytest.raises(ValueError, match="production_steps"):
        run_window(**arguments)


def test_analysis_rejects_missing_window(tmp_path: Path) -> None:
    for state in (0, 2):
        (tmp_path / f"state_{state:02d}_energies.csv").write_text(
            "sample,origin_state,u_0,u_1\n0,%d,0,1\n" % state
        )
    with pytest.raises(ValueError, match="Missing lambda windows"):
        analyze_phase(tmp_path)


def test_analysis_accepts_unequal_sample_counts(tmp_path: Path) -> None:
    for state, count in enumerate((2, 3, 5)):
        path = tmp_path / f"state_{state:02d}_energies.csv"
        with path.open("w", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(["sample", "origin_state", "u_0", "u_1", "u_2"])
            for sample in range(count):
                writer.writerow([sample, state, 0.0, 0.2, 0.4])
    result = analyze_phase(tmp_path)
    assert result["samples"] == 10
    assert result["minimum_samples_per_state"] == 2
    assert np.isfinite(result["delta_g_kj_mol"])
