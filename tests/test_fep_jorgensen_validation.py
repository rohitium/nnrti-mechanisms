from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
import pytest

from scripts.fep_jorgensen.alchemical import build_alchemical_plan, resolve_mutation_site
from scripts.fep_jorgensen.analyze import analyze_phase, analyze_target, normalize_to_reference
from scripts.fep_jorgensen.config import FEPConfig, LambdaSchedule
from scripts.fep_jorgensen.mutations import (
    MANUSCRIPT_PLANS,
    MANUSCRIPT_TARGETS,
    Mutation,
    MutationLeg,
    unique_manuscript_legs,
)
from scripts.fep_jorgensen.panel import write_worker_manifest
from scripts.fep_jorgensen.perses_hybrid import perses_available
from scripts.fep_jorgensen.worker import _perses_default_parameters, run_window


def _write_toy_phase(path: Path) -> None:
    from openmm import NonbondedForce, System, XmlSerializer, unit

    path.mkdir(parents=True)
    system = System()
    system.addParticle(39.9)
    force = NonbondedForce()
    force.addParticle(1.0 * unit.elementary_charge, 0.3 * unit.nanometer, 0.5 * unit.kilojoule_per_mole)
    system.addForce(force)
    (path / "hybrid_system.xml").write_text(XmlSerializer.serialize(system))
    (path / "hybrid_topology.pdb").write_text(
        "ATOM      1  AR  ARG A   1       0.000   0.000   0.000  1.00  0.00          Ar\n"
        "TER\nEND\n"
    )
    (path / "schedule.json").write_text(
        json.dumps(
            {
                "lambda_values": [0.0, 0.5, 1.0],
                "lambda_parameter_functions": "nonbonded-scaling",
                "alchemical_plan": {"alchemical_atom_indices": [0]},
            }
        )
    )


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


def test_resolve_v106a_site_from_md_endpoints() -> None:
    from scripts.fep_jorgensen.mutations import Mutation

    mutation = Mutation.parse("V106A")
    site = resolve_mutation_site(
        Path("results/md_runs/wt/rep_01/assets/wt_md_rep01_start.pdb"),
        Path("results/md_runs/V106A/rep_01/assets/V106A_md_rep01_start.pdb"),
        mutation,
    )
    assert site.old_residue == "VAL"
    assert site.new_residue == "ALA"
    assert site.pdb_residue_id == "103"


def test_build_alchemical_plan_for_v106a() -> None:
    from scripts.fep_jorgensen.mutations import MutationLeg

    leg = MutationLeg("WT", "V106A", "V106A")
    plan = build_alchemical_plan(leg, replicate=1)
    assert plan.strategy == "annihilate_wt_sidechain"
    assert len(plan.atom_indices) >= 6
    assert plan.start_system_xml.is_file()


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


def test_panel_prepare_commands_default_to_perses_backend() -> None:
    from scripts.fep_jorgensen.panel import preparation_commands

    command = preparation_commands(Path("results/analysis/fep_jorgensen"))[0]
    assert "--backend perses" in command


@pytest.mark.skipif(not perses_available(), reason="Perses/openmmtools not installed")
def test_openeye_shim_formats_val_template_names() -> None:
    from pkg_resources import resource_filename
    import os

    from scripts.fep_jorgensen.openeye_shim import _create_oemol_from_sdf, install_openeye_shim

    install_openeye_shim()
    pdb = resource_filename("perses", os.path.join("data", "amino_acid_templates", "VAL.pdb"))
    mol = _create_oemol_from_sdf(pdb, add_hydrogens=True)
    for atom in mol.GetAtoms():
        name = atom.GetName().replace(" ", "")
        if name and name[0].isdigit():
            atom.SetName(name[1:] + name[0])
    names = {atom.GetName() for atom in mol.GetAtoms()}
    assert {"HG11", "HG21", "CB"}.issubset(names)


@pytest.mark.skipif(not perses_available(), reason="Perses/openmmtools not installed")
def test_perses_prepare_v106a_smoke(tmp_path: Path) -> None:
    from scripts.fep_jorgensen.prepare import prepare
    from scripts.fep_jorgensen.config import FEPConfig
    from scripts.fep_jorgensen.mutations import MutationLeg

    leg = MutationLeg("WT", "V106A", "V106A")
    config = FEPConfig.for_leg(
        leg,
        output_dir=tmp_path,
        prepare_backend="perses",
    )
    prepare(config, replicate=1)
    holo = tmp_path / "legs" / leg.leg_id / "holo"
    schedule = json.loads((holo / "schedule.json").read_text())
    assert schedule["prepare_backend"] == "perses"
    assert schedule["lambda_parameter_functions"] == "perses-default"
    assert (holo / "hybrid_system.xml").is_file()


def test_panel_manifest_has_every_lambda_state_for_holo_only(tmp_path: Path) -> None:
    manifest = tmp_path / "worker_manifest.csv"
    config = FEPConfig(
        output_dir=tmp_path,
        lambda_schedule=LambdaSchedule((0.0, 0.5, 1.0)),
    )
    count = write_worker_manifest(manifest, tmp_path, config)
    rows = list(csv.DictReader(manifest.open()))
    assert count == 19 * 3
    assert len(rows) == count
    assert {row["phase"] for row in rows} == {"holo"}
    assert {int(row["state_index"]) for row in rows} == {0, 1, 2}


def test_double_mutant_analysis_sums_sequential_legs(tmp_path: Path) -> None:
    plan = MANUSCRIPT_PLANS["V106A+L234I"]
    values = ((2.0, 0.3), (4.0, 0.4))
    for leg, (delta_g, uncertainty) in zip(plan.legs, values):
        run_dir = tmp_path / "legs" / leg.leg_id
        run_dir.mkdir(parents=True)
        (run_dir / "summary.json").write_text(
            json.dumps(
                {
                    "leg_id": leg.leg_id,
                    "start_label": leg.start_label,
                    "end_label": leg.end_label,
                    "mutation": leg.mutation,
                    "delta_g_mutation_kj_mol": delta_g,
                    "uncertainty_kj_mol": uncertainty,
                }
            )
        )
    result = analyze_target("V106A_L234I", tmp_path)
    assert result["delta_g_mutation_kj_mol"] == pytest.approx(6.0)
    assert result["uncertainty_kj_mol"] == pytest.approx(0.5)


def test_reference_normalization_makes_wt_zero() -> None:
    summaries = {
        "WT": {
            "delta_g_mutation_kcal_mol": 0.5,
            "uncertainty_kcal_mol": 0.1,
        },
        "V106A": {
            "delta_g_mutation_kcal_mol": 2.5,
            "uncertainty_kcal_mol": 0.2,
        },
    }
    rows = normalize_to_reference(summaries, reference="WT")
    by_target = {row["target"]: row for row in rows}
    assert by_target["WT"]["delta_delta_g_kcal_mol"] == 0.0
    assert by_target["V106A"]["delta_delta_g_kcal_mol"] == pytest.approx(2.0)


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
    holo, windows = tmp_path / "holo", tmp_path / "holo" / "windows"
    _write_toy_phase(holo)
    for state in range(3):
        arguments = _worker_arguments(holo, windows)
        arguments["state_index"] = state
        run_window(**arguments)
    for state in range(3):
        rows = list(csv.DictReader((windows / f"state_{state:02d}_energies.csv").open()))
        assert len(rows) == 4
        assert [float(rows[0][f"u_{i}"]) for i in range(3)] == pytest.approx([0.0, 0.0, 0.0])
    result = analyze_phase(holo, temperature_k=300.0)
    assert result["samples"] == 12
    assert result["minimum_samples_per_state"] == 4
    assert result["delta_g_kj_mol"] == pytest.approx(0.0, abs=1e-6)


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
    windows = tmp_path / "windows"
    windows.mkdir()
    for state in (0, 2):
        (windows / f"state_{state:02d}_energies.csv").write_text(
            "sample,origin_state,u_0,u_1\n0,%d,0,1\n" % state
        )
    with pytest.raises(ValueError, match="Missing lambda windows"):
        analyze_phase(tmp_path)


def test_analysis_accepts_unequal_sample_counts(tmp_path: Path) -> None:
    windows = tmp_path / "windows"
    windows.mkdir()
    for state, count in enumerate((2, 3, 5)):
        path = windows / f"state_{state:02d}_energies.csv"
        with path.open("w", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(["sample", "origin_state", "u_0", "u_1", "u_2"])
            for sample in range(count):
                writer.writerow([sample, state, 0.0, 0.2, 0.4])
    result = analyze_phase(tmp_path)
    assert result["samples"] == 10
    assert result["minimum_samples_per_state"] == 2
    assert np.isfinite(result["delta_g_kj_mol"])
