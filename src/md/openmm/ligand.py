from __future__ import annotations

from pathlib import Path

from .require import require_module


def load_ligand_molecule(ligand_sdf: Path):
    openff = require_module("openff.toolkit")
    mols = openff.topology.Molecule.from_file(
        str(ligand_sdf), allow_undefined_stereo=True
    )
    molecules = mols if isinstance(mols, list) else [mols]
    if not molecules:
        raise ValueError(f"No molecules found in {ligand_sdf}")
    for mol in molecules:
        mol.assign_partial_charges(partial_charge_method="gasteiger")
    return molecules[0]


def build_forcefield(ligand_molecules) -> "openmm.app.ForceField":
    app = require_module("openmm.app")
    generators = require_module("openmmforcefields.generators")

    forcefield = app.ForceField(
        "amber14/protein.ff14SB.xml",
        "amber14/DNA.bsc1.xml",
        "amber14/tip3p.xml",
    )
    generator = generators.SMIRNOFFTemplateGenerator(
        molecules=ligand_molecules,
        template_generator_kwargs={"partial_charge_method": "gasteiger"},
    )
    forcefield.registerTemplateGenerator(generator.generator)
    return forcefield
