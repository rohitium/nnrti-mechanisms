# src/ module reference

Current code supports the OpenMM-only Sherlock workflow:
1. local asset preparation
2. cluster FEP execution
3. local result collection

## main.py
CLI entrypoint for:
- `--prepare-local-openmm-only`
- `--generate-slurm`
- `--collect-results`

## dor_alchemy_pipeline.py
Prepares WT + mutant CIFs and prebuilt OpenMM alchemical assets for cluster execution.
Key symbol: `prepare_local_openmm_only_for_cluster`.

## config.py
Run/structure specification for the active DOR/4NCG workflow.
Key symbols: `StructureSpec`, `RunSpec`, `dor_4ncg_spec`.

## cluster/
Cluster execution and postprocessing:
- `manifest.py`: task schema + CSV IO.
- `fep_worker.py`: executes one FEP leg from manifest, writes JSON + DCD outputs.
- `slurm_generator.py`: generates SLURM array script.
- `result_collector.py`: collects leg/window results, computes ΔΔG, trajectory-averaged
  structural metrics, correlations, and writes CSV outputs.

## openmm/
OpenMM integration utilities:
- `alchemy.py`: alchemical FEP engine + asset prep helpers.
- `structure.py`: restrained minimization and ligand insertion/jitter.
- `minimizer.py`: minimization utilities.
- `ligand.py`: ligand loading + forcefield template setup.
- `platform.py`: platform selection from env.
- `restraints.py`: atom selection for restraints.
- `require.py`: dependency import helper.

## analysis_metrics.py
Structural metrics from minimized structures:
contacts, H-bonds, and pocket-volume proxy.
Also includes trajectory ensemble averaging from complex-leg physical-state DCDs.

## susceptibility_io.py
Loads DOR susceptibility workbook and normalizes mutation rows.

## numbering.py
Detects auth vs label numbering used by PDBFixer.

## mutation/
Minimal mutation utilities used by preparation:
- `helpers.py`
- `mutagenesis.py`
- `steps.py`

## utils/
Shared utilities for paths, CIF parsing, residue maps, and mutation token parsing.
