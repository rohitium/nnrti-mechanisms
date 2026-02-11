"""Snakemake script: minimize structure and prepare solvated MD assets."""

import logging
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

from src.config import dor_4ncg_spec
from src.openmm.md_protocol import MDProtocolConfig, prepare_md_assets
from src.openmm.minimizer import minimize_system
from src.openmm.require import require_module
from src.openmm.structure import minimize_with_restraints

app = require_module("openmm.app")

root = Path(".").resolve()
spec = dor_4ncg_spec(root)

cif_path = Path(snakemake.input.cif)  # noqa: F821
ligand_sdf = Path(snakemake.input.ligand_sdf)  # noqa: F821
min_pdb = Path(snakemake.output.minimized_pdb)  # noqa: F821
system_xml = Path(snakemake.output.system_xml)  # noqa: F821
topology_pdb = Path(snakemake.output.topology_pdb)  # noqa: F821
seed = snakemake.params.seed  # noqa: F821
jitter = snakemake.params.jitter  # noqa: F821

# Step 1: minimize with restraints, then unrestrained polish
min_pdb.parent.mkdir(parents=True, exist_ok=True)
topology, positions, forcefield = minimize_with_restraints(
    cif_path=cif_path,
    ligand_resname=spec.structure.ligand_resname,
    ligand_sdf=ligand_sdf,
    restraint_radius_angstrom=spec.restraint_radius_angstrom,
    restraint_k_kj_mol_nm2=spec.restraint_k_kj_mol_nm2,
    output_path=min_pdb,
    jitter_seed=seed,
    jitter_angstrom=jitter,
)
_, positions = minimize_system(
    topology,
    positions,
    forcefield,
    restraint_indices=[],
    restraint_k_kj_mol_nm2=0.0,
)
with open(min_pdb, "w") as handle:
    app.PDBFile.writeFile(topology, positions, handle)

# Step 2: solvate and serialize system
system_xml.parent.mkdir(parents=True, exist_ok=True)
prepare_md_assets(
    minimized_pdb_path=min_pdb,
    ligand_resname=spec.structure.ligand_resname,
    ligand_sdf=ligand_sdf,
    topology_pdb_path=topology_pdb,
    system_xml_path=system_xml,
    config=MDProtocolConfig(),
)
