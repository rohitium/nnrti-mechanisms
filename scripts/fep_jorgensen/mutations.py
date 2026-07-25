from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re


AA1_TO_AA3 = {
    "A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS",
    "Q": "GLN", "E": "GLU", "G": "GLY", "H": "HIS", "I": "ILE",
    "L": "LEU", "K": "LYS", "M": "MET", "F": "PHE", "P": "PRO",
    "S": "SER", "T": "THR", "W": "TRP", "Y": "TYR", "V": "VAL",
}


def canonical_label(label: str) -> str:
    value = label.strip().upper().replace("_", "+")
    return "WT" if value == "WT" else value


def safe_label(label: str) -> str:
    value = canonical_label(label)
    return "wt" if value == "WT" else value.replace("+", "_")


def mutation_tokens(label: str) -> tuple[str, ...]:
    value = canonical_label(label)
    if value == "WT":
        return ()
    tokens = tuple(value.split("+"))
    for token in tokens:
        Mutation.parse(token)
    if len(set(tokens)) != len(tokens):
        raise ValueError(f"Duplicate substitutions in {label}")
    return tokens


@dataclass(frozen=True)
class Mutation:
    label: str
    residue_id: str
    old_residue: str
    new_residue: str

    @classmethod
    def parse(cls, label: str) -> "Mutation":
        value = canonical_label(label)
        match = re.fullmatch(r"([A-Z])(\d+)([A-Z])", value)
        if not match:
            raise ValueError(f"Expected one substitution such as V106A, got {label!r}")
        old, residue_id, new = match.groups()
        return cls(value, residue_id, AA1_TO_AA3[old], AA1_TO_AA3[new])


@dataclass(frozen=True)
class MutationLeg:
    start_label: str
    end_label: str
    mutation: str

    def __post_init__(self) -> None:
        object.__setattr__(self, "start_label", canonical_label(self.start_label))
        object.__setattr__(self, "end_label", canonical_label(self.end_label))
        object.__setattr__(self, "mutation", Mutation.parse(self.mutation).label)
        start = set(mutation_tokens(self.start_label))
        end = set(mutation_tokens(self.end_label))
        if end != start | {self.mutation} or self.mutation in start:
            raise ValueError(
                f"Leg {self.start_label}->{self.end_label} does not add exactly {self.mutation}"
            )

    @property
    def leg_id(self) -> str:
        return f"{safe_label(self.start_label)}_to_{safe_label(self.end_label)}"

    def input_complex_pdb(self, replicate: int = 1) -> Path:
        label = safe_label(self.start_label)
        prefix = "wt" if self.start_label == "WT" else label
        return (
            Path("results/md_runs") / label / f"rep_{replicate:02d}" / "assets"
            / f"{prefix}_md_rep{replicate:02d}_start.pdb"
        )

    def input_complex_system(self, replicate: int = 1) -> Path:
        pdb = self.input_complex_pdb(replicate)
        return pdb.with_name(pdb.name.replace("_start.pdb", "_system.xml"))

    def endpoint_complex_pdb(self, replicate: int = 1) -> Path:
        label = safe_label(self.end_label)
        return (
            Path("results/md_runs") / label / f"rep_{replicate:02d}" / "assets"
            / f"{label}_md_rep{replicate:02d}_start.pdb"
        )

    def _apo_assets_dir(self, label: str, replicate: int) -> Path:
        return (
            Path("results/md_runs") / "apo" / safe_label(label) / f"rep_{replicate:02d}" / "assets"
        )

    def input_apo_pdb(self, replicate: int = 1) -> Path:
        label = safe_label(self.start_label)
        prefix = "wt" if self.start_label == "WT" else label
        return (
            self._apo_assets_dir(self.start_label, replicate)
            / f"{prefix}_apo_md_rep{replicate:02d}_start.pdb"
        )

    def input_apo_system(self, replicate: int = 1) -> Path:
        pdb = self.input_apo_pdb(replicate)
        return pdb.with_name(pdb.name.replace("_start.pdb", "_system.xml"))

    def endpoint_apo_pdb(self, replicate: int = 1) -> Path:
        label = safe_label(self.end_label)
        return (
            self._apo_assets_dir(self.end_label, replicate)
            / f"{label}_apo_md_rep{replicate:02d}_start.pdb"
        )


@dataclass(frozen=True)
class TargetPlan:
    target: str
    legs: tuple[MutationLeg, ...]

    def __post_init__(self) -> None:
        object.__setattr__(self, "target", canonical_label(self.target))
        if not self.legs or self.legs[0].start_label != "WT":
            raise ValueError(f"Plan for {self.target} must start at WT")
        for previous, following in zip(self.legs, self.legs[1:]):
            if previous.end_label != following.start_label:
                raise ValueError(f"Discontinuous plan for {self.target}")
        if self.legs[-1].end_label != self.target:
            raise ValueError(f"Plan ends at {self.legs[-1].end_label}, not {self.target}")


def _single(label: str) -> TargetPlan:
    value = canonical_label(label)
    return TargetPlan(value, (MutationLeg("WT", value, value),))


MANUSCRIPT_PLANS: dict[str, TargetPlan] = {
    label: _single(label)
    for label in (
        "F227C", "G190A", "G190E", "G190S", "K103N", "V106A", "V106I",
        "V106M", "Y181C", "Y188L", "Y318F",
    )
}

for target, intermediate, added_mutation in (
    ("A98G+F227C", "F227C", "A98G"),
    ("K103N+M230L", "K103N", "M230L"),
    ("K103N+P225H", "K103N", "P225H"),
    ("L100I+K103N", "K103N", "L100I"),
    ("V106A+F227L", "V106A", "F227L"),
    ("V106A+L234I", "V106A", "L234I"),
    ("V106A+P225H", "V106A", "P225H"),
    ("V106I+F227C", "V106I", "F227C"),
):
    MANUSCRIPT_PLANS[target] = TargetPlan(
        target,
        (
            MANUSCRIPT_PLANS[intermediate].legs[0],
            MutationLeg(intermediate, target, added_mutation),
        ),
    )


MANUSCRIPT_TARGETS = tuple(MANUSCRIPT_PLANS)


def unique_manuscript_legs() -> tuple[MutationLeg, ...]:
    unique: dict[str, MutationLeg] = {}
    for plan in MANUSCRIPT_PLANS.values():
        for leg in plan.legs:
            unique.setdefault(leg.leg_id, leg)
    return tuple(unique.values())
