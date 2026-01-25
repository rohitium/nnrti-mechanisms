# src/ module reference

This document lists each Python module in `src/`, its purpose, inputs/outputs,
and the functions/classes defined in that module.

## analysis_metrics.py
Purpose
- Compute ligand-protein contact metrics and a pocket-volume proxy from minimized structures.

Inputs
- `pdbx_path`: path to minimized structure (`.pdb`).
- `ligand_resname`: ligand residue name.

Outputs
- Contact counts, H-bond counts, and pocket volume proxy (A^3).

Defined symbols
- `ContactMetrics` (dataclass): `contact_count`, `hbond_count`.
- `compute_contacts(pdbx_path: Path, ligand_resname: str, cutoff_angstrom: float = 4.0) -> ContactMetrics`
- `pocket_volume_proxy(pdbx_path: Path, ligand_resname: str, grid_spacing: float = 0.5, radius_angstrom: float = 8.0) -> float`

## config.py
Purpose
- Central configuration for RPV/DOR structures and default run parameters.

Inputs
- `root`: repository root path.

Outputs
- `RunSpec` objects for RPV and DOR.

Defined symbols
- `MutationSpec` (dataclass)
- `StructureSpec` (dataclass)
- `RunSpec` (dataclass)
- `rpv_spec(root: Path) -> RunSpec`
- `dor_spec(root: Path) -> RunSpec`

## ligand_from_cif.py
Purpose
- Generate hydrogenated ligand SDF files from CIF metadata.

Inputs
- CIF file, ligand residue name (comp_id), ligand chain id, output SDF path.

Outputs
- Hydrogenated ligand SDF at the requested path.

Defined symbols
- `LigandAtom` (dataclass)
- `LigandBond` (dataclass)
- `_require(module_name: str)` (internal import helper)
- `_bond_order(value_order: str) -> int`
- `_normalize_element(element: str) -> str`
- `_load_block(cif_path: Path)`
- `_chem_comp_atoms(block, comp_id: str) -> dict[str, tuple[str, bool]]`
- `_chem_comp_bonds(block, comp_id: str) -> list[LigandBond]`
- `_atom_site_category(block) -> dict`
- `_ligand_atoms(block, comp_id: str, chain_id: str, atom_elements: dict[str, tuple[str, bool]]) -> list[LigandAtom]`
- `_build_rdkit_mol(atoms: list[LigandAtom], bonds: list[LigandBond])`
- `generate_ligand_sdf(cif_path: Path, comp_id: str, chain_id: str, out_path: Path) -> Path`
- `_add_explicit_hydrogens(mol)`
- `_write_sdf_file(mol, out_path: Path) -> None`
- `main() -> None`

## main.py
Purpose
- Orchestrate the full DRM pipeline: WT minimization, mutation application,
  metrics computation, and plot generation.

Inputs
- `data/DRMs.csv` (mutation list with chain column).
- Structure/ligand paths from `config.py`.

Outputs
- `results/metrics_summary.csv`
- `results/plots/*_delta_metrics.png`

Defined symbols
- `_load_drms(drms_path: Path) -> pd.DataFrame`
- `_prepare_structure(cif_path: Path, ligand_resname: str, ligand_sdf: Path, restraint_radius: float, restraint_k: float, output_path: Path)`
- `_mutation_worker(task: dict) -> dict`
- `_run_mutations(run_spec, paths, mutation_rows: pd.DataFrame, chain_map: dict[str, str])`
- `main() -> None`

## mutagenesis.py
Purpose
- Apply one or multiple point mutations to a CIF structure using PDBFixer.

Inputs
- CIF path, chain id, residue id, new residue.

Outputs
- Mutated CIF written to output path.

Defined symbols
- `_require(module_name: str)` (internal import helper)
- `_three_letter(res_name: str) -> str`
- `_mutation_strings(old_res: str, res_id: str, new_res: str) -> Iterable[str]`
- `_residue_name_in_chain(fixer, chain_id: str, residue_id: str) -> str`
- `apply_mutation(cif_path: Path, chain_id: str, residue_id: str, new_residue: str, output_path: Path) -> Path`
- `apply_mutations(cif_path: Path, mutations: Iterable[tuple[str, str, str]], output_path: Path) -> Path`

## openmm_pipeline.py
Purpose
- Minimize complexes with OpenMM and compute energy-based binding proxy terms.

Inputs
- CIF path, ligand residue name, ligand SDF, restraint settings.

Outputs
- Minimized CIF/PDB files and energy proxy values.

Defined symbols
- `_require(module_name: str)` (internal import helper)
- `EnergyResult` (dataclass)
- `load_ligand_molecule(ligand_sdf: Path)`
- `build_forcefield(ligand_molecules) -> "openmm.app.ForceField"`
- `_minimize(topology, positions, forcefield, restraint_indices: Sequence[int], restraint_k_kj_mol_nm2: float)`
- `_modeller_without(topology, positions, residues_to_delete)`
- `_compute_energy(topology, positions, forcefield) -> float`
- `_heavy_atom_indices(topology, exclude_resname: str) -> list[int]`
- `_restrained_indices(positions, ligand_indices: Sequence[int], candidate_indices: Sequence[int], radius_angstrom: float) -> list[int]`
- `minimize_with_restraints(cif_path: Path, ligand_resname: str, ligand_sdf: Path, restraint_radius_angstrom: float, restraint_k_kj_mol_nm2: float, output_path: Path) -> tuple`
- `_get_platform()`
- `compute_binding_proxy(topology, positions, forcefield, ligand_resname: str) -> EnergyResult`

## plotting.py
Purpose
- Generate per-drug delta metric bar charts.

Inputs
- DataFrame with metrics and a `Paths` object for output locations.

Outputs
- `results/plots/*_delta_metrics.png`

Defined symbols
- `plot_delta_metrics(df: pd.DataFrame, paths) -> None`

## utils.py
Purpose
- Shared helpers for filesystem paths, label parsing, and CIF chain mapping.

Inputs/Outputs
- See function signatures below.

Defined symbols
- `Paths` (dataclass)
- `project_paths(root: Path) -> Paths`
- `ensure_dirs(paths: Iterable[Path]) -> None`
- `sanitize_label(label: str) -> str`
- `parse_mutation_token(token: str) -> tuple[str, str, str]`
- `parse_mutation_group(mutation: str, chains: str | list[str]) -> list[tuple[str, str, str]]`
- `load_chain_subunits(cif_path: Path) -> dict[str, str]`
