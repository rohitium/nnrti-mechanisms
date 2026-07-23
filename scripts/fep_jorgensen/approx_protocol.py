"""OpenMM/openmmtools approximation of the Jorgensen (2000) mutation FEP cycle.

This is intentionally *not* an exact MCPRO reproduction.  It keeps the paper's
scientific shape — MD equilibration of each inhibitor complex, alchemical
protein-side-chain mutation in the bound complex, MBAR estimation, and
inhibitor-relative normalization — while using the manuscript's OpenMM stack.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
import json
from pathlib import Path


@dataclass(frozen=True)
class ApproxJorgensenProtocol:
    """Machine-readable contract for the merged approximate workflow."""

    label: str = (
        "OpenMM/openmmtools MD-equilibrated MCMC alchemical mutation FEP, "
        "Jorgensen-cycle inspired"
    )
    force_field_protein: str = "amber14/protein.ff14SB.xml"
    force_field_dna: str = "amber14/DNA.bsc1.xml"
    force_field_water: str = "amber14/tip3p.xml"
    small_molecule_forcefield: str = "openff-2.0.0"
    solvent_model: str = "TIP3P explicit PME"
    md_integrator: str = "LangevinMiddleIntegrator"
    md_timestep_fs: float = 2.0
    md_minimization_iterations: int = 500
    md_initial_temperature_k: float = 100.0
    md_initial_equilibration_ps: float = 3.0
    md_equilibration_temperature_k: float = 300.0
    md_equilibration_ps: float = 50.0
    quench_blocks: int = 6
    quench_block_ps: float = 4.0
    quench_start_k: float = 300.0
    quench_end_k: float = 50.0
    mc_sampler: str = "fixed-lambda Langevin windows with multistate energy reevaluation"
    mc_alternative_sampler: str = "openmmtools GHMC/HMC multistate replica exchange"
    lambda_windows: int = 11
    equilibration_steps_per_window: int = 250_000
    production_steps_per_window: int = 2_500_000
    energy_interval_steps: int = 2_500
    thermodynamic_cycle: str = "protein-side-chain mutation in inhibitor-bound complex only"
    reference_system: str = "WT"
    manuscript_inhibitor: str = "doravirine"
    exact_protocol_module: str = "scripts.fep_jorgensen.exact_protocol"

    def validate(self) -> None:
        if self.md_timestep_fs <= 0:
            raise ValueError("md_timestep_fs must be positive")
        if self.lambda_windows < 3:
            raise ValueError("lambda_windows must be at least 3")
        if self.production_steps_per_window % self.energy_interval_steps:
            raise ValueError("production_steps_per_window must divide energy_interval_steps")

    @property
    def md_initial_equilibration_steps(self) -> int:
        return self._ps_to_steps(self.md_initial_equilibration_ps)

    @property
    def md_equilibration_steps(self) -> int:
        return self._ps_to_steps(self.md_equilibration_ps)

    @property
    def quench_block_steps(self) -> int:
        return self._ps_to_steps(self.quench_block_ps)

    def _ps_to_steps(self, duration_ps: float) -> int:
        steps = int(round(duration_ps * 1000.0 / self.md_timestep_fs))
        if steps < 1:
            raise ValueError(f"Duration {duration_ps} ps is too short for timestep")
        return steps

    def to_dict(self) -> dict:
        self.validate()
        payload = asdict(self)
        payload.update(
            {
                "md_initial_equilibration_steps": self.md_initial_equilibration_steps,
                "md_equilibration_steps": self.md_equilibration_steps,
                "quench_block_steps": self.quench_block_steps,
                "lambda_values": [i / (self.lambda_windows - 1) for i in range(self.lambda_windows)],
            }
        )
        return payload

    def write(self, path: Path) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(self.to_dict(), indent=2) + "\n")
