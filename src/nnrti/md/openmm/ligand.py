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


#: SMIRNOFF release used to parameterise the ligand; verified against the
#: production system XMLs. See build_forcefield for why this is pinned.
SMIRNOFF_FORCEFIELD = "openff-2.0.0"


def build_forcefield(ligand_molecules) -> "openmm.app.ForceField":
    app = require_module("openmm.app")
    generators = require_module("openmmforcefields.generators")

    forcefield = app.ForceField(
        "amber14/protein.ff14SB.xml",
        "amber14/DNA.bsc1.xml",
        "amber14/tip3p.xml",
    )
    # Pinned explicitly. openmmforcefields <0.16 resolved forcefield=None to a
    # default release; 0.16 removed that fallback and raises instead, so the
    # version has to be named. openff-2.0.0 reproduces the ligand parameters
    # recorded in the production system XMLs bit-identically (sigma, epsilon and
    # charge all match to 0.0); openff-2.0.0 through 2.2.1 are equivalent for
    # doravirine, while 2.3.0 differs substantially. Changing this silently
    # reparameterises the ligand away from the MD that produced the trajectories.
    generator = generators.SMIRNOFFTemplateGenerator(
        molecules=ligand_molecules,
        forcefield=SMIRNOFF_FORCEFIELD,
        template_generator_kwargs={"partial_charge_method": "gasteiger"},
    )
    forcefield.registerTemplateGenerator(generator.generator)
    return forcefield
