"""Machine-readable specification of the Rizzo/Jorgensen (2000) FEP protocol.

This module deliberately contains no OpenMM/Perses substitutions.  An execution
backend may call IMPACT c1.00 and MCPRO 1.65, but it must satisfy every invariant
validated here before its results are labelled an exact protocol reproduction.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
import json
from pathlib import Path
from typing import Mapping


MODEL_RESIDUES = (
    "91-110A", "161-205A", "222-242A", "316-321A", "343-349A",
    "381-383A", "134-140B",
)
MD_FREE_RESIDUES = (
    "95-107A", "172A", "177-182A", "188-192A", "198A", "227A",
    "229A", "234-236A", "318-319A", "321A", "135-139B",
)
MD_RESTRAINED_RESIDUES = (
    "94A", "108A", "175-176A", "183A", "187A", "225A",
    "237-239A", "317A", "320A", "349A", "382-383A", "134B", "140B",
)
MC_RIGID_RESIDUES = (
    "91-94A", "109-110A", "116-178A", "184-185A", "192-197A",
    "199-205A", "222-224A", "230-232A", "240-242A", "316-317A",
    "320-321A", "343-349A", "381-383A", "134-135B", "137B", "140B",
)
MC_FLEXIBLE_RESIDUES = (
    "95-108A", "179-183A", "186-191A", "198A", "225-229A",
    "233-239A", "318-319A", "136B", "138B",
)
SOURCE_STRUCTURES: Mapping[str, tuple[str, ...]] = {
    "sustiva": ("1rt1",),
    "nevirapine": ("1vrt",),
    "mkc-442": ("1rti", "1rt1"),
    "9-cl-tibo": ("1rev",),
}


@dataclass(frozen=True)
class ExactJorgensenProtocol:
    impact_version: str = "c1.00"
    mcpro_version: str = "1.65"
    force_field: str = "CM1P-augmented OPLS-AA"
    dielectric: str = "epsilon=4r"
    md_minimization_cycles: int = 10
    md_integrator: str = "Verlet"
    md_timestep_ps: float = 0.001
    thermostat: str = "Berendsen"
    thermostat_coupling_ps: float = 0.2
    constraints: str = "SHAKE bond lengths"
    md_initial_temperature_k: float = 100.0
    md_initial_equilibration_ps: float = 3.0
    md_equilibration_temperature_k: float = 300.0
    md_equilibration_ps: float = 50.0
    quench_blocks: int = 6
    quench_block_ps: float = 4.0
    quench_start_k: float = 300.0
    quench_end_k: float = 50.0
    water_model: str = "TIP4P"
    water_cap_radius_angstrom: float = 22.0
    expected_water_count: int = 850
    mc_solvent_equilibration_configurations: int = 1_000_000
    mc_full_equilibration_configurations: int = 10_000_000
    mc_averaging_configurations: int = 10_000_000
    backbone_fixed_during_mc: bool = True
    inhibitor_fully_flexible_during_mc: bool = True
    reference_inhibitor: str = "sustiva"

    def validate(self) -> None:
        required = {
            "impact_version": "c1.00",
            "mcpro_version": "1.65",
            "force_field": "CM1P-augmented OPLS-AA",
            "dielectric": "epsilon=4r",
            "md_integrator": "Verlet",
            "md_timestep_ps": 0.001,
            "thermostat": "Berendsen",
            "thermostat_coupling_ps": 0.2,
            "constraints": "SHAKE bond lengths",
            "water_model": "TIP4P",
            "water_cap_radius_angstrom": 22.0,
            "mc_solvent_equilibration_configurations": 1_000_000,
            "mc_full_equilibration_configurations": 10_000_000,
            "mc_averaging_configurations": 10_000_000,
            "backbone_fixed_during_mc": True,
            "inhibitor_fully_flexible_during_mc": True,
            "reference_inhibitor": "sustiva",
        }
        mismatches = [
            f"{name}={getattr(self, name)!r} (required {value!r})"
            for name, value in required.items() if getattr(self, name) != value
        ]
        if (self.md_initial_equilibration_ps, self.md_equilibration_ps) != (3.0, 50.0):
            mismatches.append("MD equilibration must be 3 ps at 100 K then 50 ps at 300 K")
        if (self.quench_blocks, self.quench_block_ps, self.quench_start_k, self.quench_end_k) != (6, 4.0, 300.0, 50.0):
            mismatches.append("quench must be six 4 ps blocks from 300 K to 50 K")
        if mismatches:
            raise ValueError("Not the exact Jorgensen protocol: " + "; ".join(mismatches))

    def to_dict(self) -> dict:
        self.validate()
        return {
            **asdict(self),
            "model_residues": MODEL_RESIDUES,
            "md_free_residues": MD_FREE_RESIDUES,
            "md_restrained_residues": MD_RESTRAINED_RESIDUES,
            "mc_rigid_residues": MC_RIGID_RESIDUES,
            "mc_flexible_residues": MC_FLEXIBLE_RESIDUES,
            "source_structures": SOURCE_STRUCTURES,
        }

    def write(self, path: Path) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(self.to_dict(), indent=2) + "\n")
