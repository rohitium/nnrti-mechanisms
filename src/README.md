# src/ module reference

Each module in `src/` is listed with purpose, inputs/outputs, and key symbols.

## __init__.py
Purpose: package marker for the NNRTI mutagenesis pipeline.
Symbols: (none; module docstring only).

## analysis_metrics.py
Purpose: compute contacts, H-bonds, and pocket-volume proxy.
Inputs: minimized `.pdb`, ligand resname. Outputs: `ContactMetrics`, volume.
Symbols: `ContactMetrics`, `compute_contacts`, `pocket_volume_proxy`.

## config.py
Purpose: structure/run specs for RPV/DOR.
Symbols: `MutationSpec`, `StructureSpec`, `RunSpec`, `rpv_spec`, `dor_spec`.

## drm_io.py
Purpose: load and normalize `data/DRMs.csv`.
Symbols: `load_drms`.

## main.py
Purpose: CLI entrypoint for validation and pipeline runs.
Inputs: DRM CSV, structures, ligands. Outputs: metrics + plots.
Symbols: `main`.

## metrics_io.py
Purpose: write `metrics_summary.xlsx` (RPV/DOR sheets).
Symbols: `write_metrics_xlsx`.

## numbering.py
Purpose: detect auth vs label numbering used by PDBFixer.
Symbols: `detect_numbering_scheme`.

## plotting.py
Purpose: generate per-drug delta bar plots.
Symbols: `plot_delta_metrics`.

## structure_prep.py
Purpose: restrained minimization, unrestrained minimization, and metrics.
Symbols: `prepare_structure`.

## validation.py
Purpose: validate DRM substitutions and optionally verify mutations without OpenMM.
Symbols: `validate_mutations`, `verify_mutations_only`.

## ligand_cif/
Purpose: CIF ligand parsing and SDF generation.
- `ligand_cif/types.py`: data classes for CIF ligand parsing (`LigandAtom`, `LigandBond`).
- `ligand_cif/block.py`: CIF block utilities (`bond_order`, `normalize_element`, `load_block`).
- `ligand_cif/comp.py`: parse chem_comp atoms/bonds (`chem_comp_atoms`, `chem_comp_bonds`).
- `ligand_cif/atoms.py`: parse ligand atom records from `_atom_site` (`atom_site_category`, `ligand_atoms`).
- `ligand_cif/build.py`: build RDKit molecule and write SDF (`build_rdkit_mol`, `add_explicit_hydrogens`, `write_sdf_file`).
- `ligand_cif/from_cif.py`: CLI to generate ligand SDFs (`generate_ligand_sdf`, `main`).

## mutation/
Purpose: mutation parsing, application, verification, and orchestration.
- `mutation/helpers.py`: shared helpers for PDBFixer mutations.
- `mutation/mutagenesis.py`: apply one or multiple mutations to CIF.
- `mutation/steps.py`: validate DRM tokens and map to structure residue IDs.
- `mutation/tasks.py`: build mutation tasks (including subset expansions) and compute WT metrics.
- `mutation/rows.py`: assemble metrics rows from results.
- `mutation/runner.py`: run mutation tasks in parallel and collect results.
- `mutation/worker.py`: apply a mutation task and compute metrics.
- `mutation/verify.py`: verify applied mutations by comparing base and mutated CIFs.

## openmm/
Purpose: OpenMM integration utilities.
- `openmm/require.py`: shared import helper for OpenMM stack.
- `openmm/platform.py`: select OpenMM platform and properties.
- `openmm/ligand.py`: load ligand and build forcefield templates.
- `openmm/restraints.py`: atom selection for restraints.
- `openmm/minimizer.py`: run OpenMM minimization utilities (implicit MD helper optional).
- `openmm/structure.py`: prepare complex + minimize with restraints.
- `openmm/energy.py`: compute binding proxy energies.
- `openmm/pipeline.py`: public OpenMM API re-exports.

## utils/
Purpose: shared utilities for CIF parsing, mutations, and paths.
- `utils/cif_parser.py`: mmCIF loop tokenizer for local parsers.
- `utils/cif.py`: chain→subunit mapping from CIF.
- `utils/residue_map.py`: build auth/label residue maps from CIF atom records.
- `utils/mutations.py`: mutation token parsing and label normalization.
- `utils/paths.py`: project paths and directory creation.
- `utils/__init__.py`: convenience re-exports (`Paths`, `project_paths`, `ensure_dirs`,
  `load_chain_subunits`, `load_residue_mappings`, `sanitize_label`,
  `parse_mutation_token`, `one_to_three`).
