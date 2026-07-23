from __future__ import annotations

import csv
from pathlib import Path

import pytest

from scripts.fep_jorgensen.analyze_exact import normalize_to_sustiva, read_mutation_legs
from scripts.fep_jorgensen.exact_protocol import (
    ExactJorgensenProtocol,
    MC_FLEXIBLE_RESIDUES,
    MC_RIGID_RESIDUES,
    MD_FREE_RESIDUES,
    MD_RESTRAINED_RESIDUES,
    MODEL_RESIDUES,
)


def test_exact_protocol_matches_published_constants() -> None:
    protocol = ExactJorgensenProtocol()
    protocol.validate()
    data = protocol.to_dict()
    assert data["impact_version"] == "c1.00"
    assert data["mcpro_version"] == "1.65"
    assert data["force_field"] == "CM1P-augmented OPLS-AA"
    assert data["dielectric"] == "epsilon=4r"
    assert data["md_timestep_ps"] == 0.001
    assert data["water_model"] == "TIP4P"
    assert data["water_cap_radius_angstrom"] == 22.0
    assert data["model_residues"] == MODEL_RESIDUES
    assert data["md_free_residues"] == MD_FREE_RESIDUES
    assert data["md_restrained_residues"] == MD_RESTRAINED_RESIDUES
    assert data["mc_rigid_residues"] == MC_RIGID_RESIDUES
    assert data["mc_flexible_residues"] == MC_FLEXIBLE_RESIDUES


def test_exact_protocol_rejects_openmm_style_substitution() -> None:
    with pytest.raises(ValueError, match="Not the exact Jorgensen protocol"):
        ExactJorgensenProtocol(water_model="TIP3P").validate()
    with pytest.raises(ValueError, match="Not the exact Jorgensen protocol"):
        ExactJorgensenProtocol(md_integrator="LangevinMiddle").validate()


def test_jorgensen_cycle_is_normalized_to_sustiva() -> None:
    rows = normalize_to_sustiva(
        {
            "sustiva": (1.25, 0.20),
            "nevirapine": (3.00, 0.30),
        },
        "V106A",
    )
    by_inhibitor = {row["inhibitor"]: row for row in rows}
    assert by_inhibitor["sustiva"]["delta_delta_g_kcal_mol"] == 0.0
    assert by_inhibitor["sustiva"]["uncertainty_kcal_mol"] == 0.0
    assert by_inhibitor["nevirapine"]["delta_delta_g_kcal_mol"] == pytest.approx(1.75)
    assert by_inhibitor["nevirapine"]["uncertainty_kcal_mol"] == pytest.approx(
        (0.20**2 + 0.30**2) ** 0.5
    )


def test_exact_analysis_requires_sustiva(tmp_path: Path) -> None:
    path = tmp_path / "legs.csv"
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["inhibitor", "delta_g_mutation_kcal_mol", "uncertainty_kcal_mol"])
        writer.writerow(["nevirapine", 2.0, 0.3])
    with pytest.raises(ValueError, match="requires a sustiva leg"):
        read_mutation_legs(path)
