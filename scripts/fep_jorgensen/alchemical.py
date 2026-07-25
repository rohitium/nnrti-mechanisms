"""Mutation-agnostic alchemical site setup from existing MD assets.

We already have equilibrated holo complexes and serialized OpenMM systems for
every genotype.  FEP does not need new MD or Perses hybrids: it needs to know
which atoms in the *start* system should be scaled alchemically along a path
toward the *end* genotype.

Site mapping uses endpoint PDB comparison (ordinal in the protein chain).
Atom selection uses Amber ff14SB residue templates from OpenMM.
"""

from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path

from .mutations import Mutation, MutationLeg

PROTEIN_RESIDUES = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "HID",
    "HIE", "HIP", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR",
    "VAL", "TRP", "TYR",
}
BACKBONE_ATOMS = {"N", "H", "CA", "HA", "HA2", "HA3", "C", "O"}
HIS_VARIANTS = {"HIS", "HID", "HIE", "HIP"}


@dataclass(frozen=True)
class MutationSite:
    chain_id: str
    ordinal: int
    pdb_residue_id: str
    old_residue: str
    new_residue: str
    mutation: str


@dataclass(frozen=True)
class AlchemicalPlan:
    site: MutationSite
    strategy: str
    atom_indices: tuple[int, ...]
    atom_names: tuple[str, ...]
    start_pdb: Path
    start_system_xml: Path
    end_pdb: Path
    annihilate_count: int
    insert_count: int

    def to_dict(self) -> dict:
        return {
            "strategy": self.strategy,
            "chain_id": self.site.chain_id,
            "ordinal": self.site.ordinal,
            "pdb_residue_id": self.site.pdb_residue_id,
            "old_residue": self.site.old_residue,
            "new_residue": self.site.new_residue,
            "mutation": self.site.mutation,
            "alchemical_atom_indices": list(self.atom_indices),
            "alchemical_atom_names": list(self.atom_names),
            "start_pdb": str(self.start_pdb),
            "start_system_xml": str(self.start_system_xml),
            "end_pdb": str(self.end_pdb),
            "annihilate_template_atoms": self.annihilate_count,
            "insert_template_atoms": self.insert_count,
        }


def _protein_residues(topology, chain_id: str) -> list:
    for chain in topology.chains():
        if chain.id != chain_id:
            continue
        return [res for res in chain.residues() if res.name in PROTEIN_RESIDUES]
    raise ValueError(f"Chain {chain_id!r} not found in topology.")


def _residue_templates() -> dict[str, set[str]]:
    from openmm import app

    ff = app.ForceField("amber14/protein.ff14SB.xml")
    return {name: {atom.name for atom in template.atoms} for name, template in ff._templates.items()}


def _template_name(residue_name: str) -> str:
    if residue_name in HIS_VARIANTS:
        return residue_name if residue_name in {"HID", "HIE", "HIP"} else "HIP"
    return residue_name


def resolve_mutation_site(
    start_pdb: Path,
    end_pdb: Path,
    mutation: Mutation,
    chain_id: str = "A",
) -> MutationSite:
    """Locate the mutated residue by diffing prepared holo endpoint PDBs."""
    from openmm import app

    start_top = app.PDBFile(str(start_pdb)).topology
    end_top = app.PDBFile(str(end_pdb)).topology
    start_residues = _protein_residues(start_top, chain_id)
    end_residues = _protein_residues(end_top, chain_id)
    if len(start_residues) != len(end_residues):
        raise ValueError(
            f"Protein residue count mismatch between {start_pdb} and {end_pdb}: "
            f"{len(start_residues)} vs {len(end_residues)}"
        )
    matches = []
    for ordinal, (start_res, end_res) in enumerate(zip(start_residues, end_residues)):
        if start_res.name == mutation.old_residue and end_res.name == mutation.new_residue:
            matches.append((ordinal, start_res, end_res))
    if len(matches) != 1:
        raise ValueError(
            f"Expected exactly one {mutation.old_residue}->{mutation.new_residue} site "
            f"for {mutation.label} between {start_pdb.name} and {end_pdb.name}; "
            f"found {len(matches)}"
        )
    ordinal, start_res, end_res = matches[0]
    return MutationSite(
        chain_id=chain_id,
        ordinal=ordinal,
        pdb_residue_id=str(start_res.id),
        old_residue=start_res.name,
        new_residue=end_res.name,
        mutation=mutation.label,
    )


def _atom_name_set(residue) -> set[str]:
    return {atom.name for atom in residue.atoms()}


def _indices_for_names(residue, names: set[str]) -> tuple[int, ...]:
    atoms = {atom.name: atom.index for atom in residue.atoms()}
    missing = sorted(names.difference(atoms))
    if missing:
        raise ValueError(f"Residue {residue} lacks expected atoms: {missing}")
    return tuple(sorted(atoms[name] for name in names))


def build_alchemical_plan(
    leg: MutationLeg,
    replicate: int = 1,
    chain_id: str = "A",
) -> AlchemicalPlan:
    """Build an alchemical atom list from MD assets for one leg."""
    from openmm import app

    mutation = Mutation.parse(leg.mutation)
    start_pdb = leg.input_complex_pdb(replicate)
    start_system = leg.input_complex_system(replicate)
    end_pdb = leg.endpoint_complex_pdb(replicate)
    if not start_pdb.is_file():
        raise FileNotFoundError(start_pdb)
    if not start_system.is_file():
        raise FileNotFoundError(start_system)
    if not end_pdb.is_file():
        raise FileNotFoundError(end_pdb)

    site = resolve_mutation_site(start_pdb, end_pdb, mutation, chain_id=chain_id)
    templates = _residue_templates()
    old_template = templates[_template_name(site.old_residue)]
    new_template = templates[_template_name(site.new_residue)]
    annihilate_names = (old_template - new_template) - BACKBONE_ATOMS
    insert_names = (new_template - old_template) - BACKBONE_ATOMS

    start_top = app.PDBFile(str(start_pdb)).topology
    end_top = app.PDBFile(str(end_pdb)).topology
    start_residue = _protein_residues(start_top, chain_id)[site.ordinal]
    end_residue = _protein_residues(end_top, chain_id)[site.ordinal]

    if annihilate_names:
        strategy = "annihilate_wt_sidechain"
        use_pdb = start_pdb
        use_system = start_system
        use_residue = start_residue
        selected_names = annihilate_names.intersection(_atom_name_set(use_residue))
        if not selected_names:
            raise ValueError(
                f"No annihilation atoms from template diff found on {use_residue} "
                f"for {leg.mutation}"
            )
    elif insert_names:
        strategy = "annihilate_mutant_sidechain"
        use_pdb = end_pdb
        use_system = end_pdb.with_name(
            end_pdb.name.replace("_start.pdb", "_system.xml")
        )
        use_residue = end_residue
        selected_names = insert_names.intersection(_atom_name_set(use_residue))
        if not selected_names:
            raise ValueError(
                f"No decouplable side-chain atoms from template diff found on {use_residue} "
                f"for {leg.mutation}"
            )
    else:
        shared_sidechain = (
            (old_template & new_template)
            - BACKBONE_ATOMS
            - {"CB", "HB", "HB1", "HB2", "HB3"}
        )
        strategy = "annihilate_shared_sidechain"
        use_pdb = start_pdb
        use_system = start_system
        use_residue = start_residue
        selected_names = shared_sidechain.intersection(_atom_name_set(use_residue))
        if not selected_names:
            raise ValueError(f"No alchemical atoms identified for {leg.mutation}")

    atom_indices = _indices_for_names(use_residue, selected_names)
    return AlchemicalPlan(
        site=site,
        strategy=strategy,
        atom_indices=atom_indices,
        atom_names=tuple(sorted(selected_names)),
        start_pdb=use_pdb,
        start_system_xml=use_system,
        end_pdb=end_pdb,
        annihilate_count=len(annihilate_names),
        insert_count=len(insert_names),
    )


def write_holo_phase(
    plan: AlchemicalPlan,
    holo_dir: Path,
    lambda_values: tuple[float, ...],
    leg: MutationLeg,
) -> None:
    """Serialize OpenMM inputs for the fixed-lambda worker."""
    import shutil

    from openmm import XmlSerializer, app

    holo_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy2(plan.start_system_xml, holo_dir / "hybrid_system.xml")
    shutil.copy2(plan.start_pdb, holo_dir / "hybrid_topology.pdb")
    schedule = {
        "phase": "holo",
        "mutation": leg.mutation,
        "start_label": leg.start_label,
        "end_label": leg.end_label,
        "leg_id": leg.leg_id,
        "lambda_values": list(lambda_values),
        "lambda_parameter_functions": "nonbonded-scaling",
        "thermodynamic_cycle": "protein-side-chain mutation in inhibitor-bound complex",
        "alchemical_plan": plan.to_dict(),
    }
    (holo_dir / "schedule.json").write_text(json.dumps(schedule, indent=2) + "\n")
    (holo_dir / "alchemical_plan.json").write_text(
        json.dumps(plan.to_dict(), indent=2) + "\n"
    )
