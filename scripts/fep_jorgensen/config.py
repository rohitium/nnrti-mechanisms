from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
import json


@dataclass(frozen=True)
class LambdaSchedule:
    values: tuple[float, ...] = tuple(i / 10 for i in range(11))

    def validate(self) -> None:
        if len(self.values) < 3 or self.values[0] != 0.0 or self.values[-1] != 1.0:
            raise ValueError("Lambda schedule must span 0.0 to 1.0 with at least 3 states")
        if any(b <= a for a, b in zip(self.values, self.values[1:])):
            raise ValueError("Lambda values must be strictly increasing")


@dataclass(frozen=True)
class FEPConfig:
    """Configuration for WT(λ=0)→V106A(λ=1) in holo and apo RT."""

    mutation: str = "V106A"
    chain_id: str = "A"
    residue_id: str = "106"
    wt_residue: str = "VAL"
    mutant_residue: str = "ALA"
    wt_complex_pdb: Path = Path(
        "results/md_runs/wt/rep_01/assets/wt_md_rep01_start.pdb"
    )
    ligand_sdf: Path = Path("data/ligands/dor.sdf")
    ligand_resname: str = "2KW"
    output_dir: Path = Path("results/analysis/fep_jorgensen")
    temperature_k: float = 300.0
    pressure_atm: float = 1.0
    lambda_schedule: LambdaSchedule = LambdaSchedule()
    equilibration_steps: int = 250_000
    production_steps: int = 2_500_000
    energy_interval: int = 2_500
    checkpoint_interval: int = 25_000
    timestep_fs: float = 2.0
    collision_rate_per_ps: float = 1.0
    platform: str = "CUDA"

    @property
    def run_dir(self) -> Path:
        return self.output_dir / self.mutation

    def validate(self, require_inputs: bool = False) -> None:
        if (self.mutation, self.wt_residue, self.mutant_residue) != (
            "V106A", "VAL", "ALA"
        ):
            raise ValueError("Initial implementation is restricted to WT→V106A")
        self.lambda_schedule.validate()
        if min(self.equilibration_steps, self.production_steps, self.energy_interval) < 1:
            raise ValueError("Step counts and energy interval must be positive")
        if self.production_steps % self.energy_interval:
            raise ValueError("production_steps must be divisible by energy_interval")
        if require_inputs:
            for path in (self.wt_complex_pdb, self.ligand_sdf):
                if not path.is_file():
                    raise FileNotFoundError(path)

    def write(self, path: Path) -> None:
        data = asdict(self)
        for key in ("wt_complex_pdb", "ligand_sdf", "output_dir"):
            data[key] = str(data[key])
        data["lambda_schedule"] = list(self.lambda_schedule.values)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(data, indent=2) + "\n")
