from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
import json

from .approx_protocol import ApproxJorgensenProtocol
from .mutations import Mutation, MutationLeg


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
    """Configuration for one single-residue alchemical leg."""

    mutation: str = "V106A"
    start_label: str = "WT"
    end_label: str = "V106A"
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
    approx_protocol: ApproxJorgensenProtocol = ApproxJorgensenProtocol()
    skip_equilibration: bool = True
    prepare_backend: str = "perses"

    @classmethod
    def for_leg(cls, leg: MutationLeg, **overrides) -> "FEPConfig":
        mutation = Mutation.parse(leg.mutation)
        return cls(
            mutation=mutation.label,
            start_label=leg.start_label,
            end_label=leg.end_label,
            residue_id=mutation.residue_id,
            wt_residue=mutation.old_residue,
            mutant_residue=mutation.new_residue,
            **overrides,
        )

    @property
    def leg(self) -> MutationLeg:
        return MutationLeg(self.start_label, self.end_label, self.mutation)

    @property
    def run_dir(self) -> Path:
        return self.output_dir / "legs" / self.leg.leg_id

    @property
    def equilibrated_complex_pdb(self) -> Path:
        return self.run_dir / "inputs" / "equilibrated_complex.pdb"

    @property
    def preparation_complex_pdb(self) -> Path:
        if self.skip_equilibration:
            return self.wt_complex_pdb
        equilibrated = self.equilibrated_complex_pdb
        return equilibrated if equilibrated.is_file() else self.wt_complex_pdb

    def validate(self, require_inputs: bool = False) -> None:
        parsed = Mutation.parse(self.mutation)
        if (self.residue_id, self.wt_residue, self.mutant_residue) != (
            parsed.residue_id, parsed.old_residue, parsed.new_residue
        ):
            raise ValueError("Mutation label and residue fields are inconsistent")
        self.lambda_schedule.validate()
        if min(self.equilibration_steps, self.production_steps, self.energy_interval) < 1:
            raise ValueError("Step counts and energy interval must be positive")
        if self.production_steps % self.energy_interval:
            raise ValueError("production_steps must be divisible by energy_interval")
        if self.checkpoint_interval < 1 or self.checkpoint_interval % self.energy_interval:
            raise ValueError("checkpoint_interval must be a positive multiple of energy_interval")
        if require_inputs:
            for path in (self.wt_complex_pdb, self.ligand_sdf):
                if not path.is_file():
                    raise FileNotFoundError(path)

    def write(self, path: Path) -> None:
        data = asdict(self)
        for key in ("wt_complex_pdb", "ligand_sdf", "output_dir"):
            data[key] = str(data[key])
        data["lambda_schedule"] = list(self.lambda_schedule.values)
        data["leg_id"] = self.leg.leg_id
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(data, indent=2) + "\n")
