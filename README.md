# NNRTI Resistance Mechanisms

Understanding structural mechanisms of Non-Nucleoside Reverse Transcriptase Inhibitor (NNRTI) resistance in HIV-1.

This repository contains code, data, and analysis for evaluating the structural
impact of clinically relevant NNRTI resistance mutations on Rilpivirine (RPV)
and Doravirine (DOR) using publicly available cryo-EM RT/DNA/drug structures.

## Workflows

### 1. Local Analysis Pipeline (CPU)

For quick structural analysis on a local machine:

```bash
uv run python -m src.main
```

This runs the DRM batch pipeline using `data/DRMs.csv`, writes `results/metrics_summary.csv`,
and generates per-drug plots in `results/plots/`.

### 2. Cluster FEP Pipeline (GPU)

For rigorous alchemical free energy perturbation (FEP) calculations on Stanford Sherlock:

```
LOCAL PHASE (CPU):
  CIF → Mutagenesis → Minimization → Contacts/H-bonds/Pocket → Manifest
        (PDBFixer)    (OpenMM CPU)   (MDAnalysis)              (CSV)

CLUSTER PHASE (GPU):
  Manifest → SLURM Array Job → FEP per (mutation, replicate, leg) → JSON

AGGREGATION PHASE (CPU):
  JSON results → Collect → ΔΔG = ΔG_mut - ΔG_wt → Correlation Analysis
```

#### Step 1: Local Preparation
```bash
uv run python -m src.main --prepare-local \
    --replicates 3 --seed 42 --jitter-angstrom 0.1
```

Outputs:
- `data/prepared/dor_4ncg/*.pdb` - Minimized structures
- `results/structural_metrics.csv` - Contact counts, H-bonds, pocket volumes
- `results/fep_manifest.csv` - Task assignments (84 tasks for 14 structures × 3 replicates × 2 legs)

#### Step 2: Generate SLURM Script
```bash
uv run python -m src.main --generate-slurm
```

Output: `scripts/sherlock/submit_fep.sh`

#### Step 3: Transfer to Sherlock
```bash
rsync -avz --exclude='.venv' --exclude='.git' \
    . sherlock:/scratch/users/$USER/nnrti-mechanisms/
```

#### Step 4: Submit on Sherlock
```bash
sbatch scripts/sherlock/submit_fep.sh
```

#### Step 5: Transfer Results Back
```bash
rsync -avz sherlock:/scratch/users/$USER/nnrti-mechanisms/results/fep_runs/ \
    results/fep_runs/
```

#### Step 6: Collect and Analyze
```bash
uv run python -m src.main --collect-results
```

Outputs:
- `results/ddg_summary.csv` - ΔΔG per mutation
- `results/correlation_analysis.csv` - Pearson/Spearman correlations with fold-reduction
- `results/plots/ddg_vs_fold_reduction.png`

## Key Inputs

### Structures
- `data/structures/7Z2D.cif` - RT/DNA/RPV complex
- `data/structures/7Z2G.cif` - RT/DNA/DOR complex
- `data/structures/4NCG.cif` - RT/DOR complex (for FEP workflow)

### DRM Data
- `data/DRMs.csv` - Drug resistance mutations (from Stanford HIVDB)
- `data/DRM-susceptibilities.csv.xlsx` - DOR susceptibility data with fold-reduction values

### DRM Table Format

`data/DRMs.csv` must include:
- `drug`: `RPV` or `DOR`
- `mutation`: one or more mutations, e.g. `K101E` or `K101E+E138K`
- `chain`: chain id(s) for each mutation, e.g. `A` or `A+B`
- `category`, `notes`: optional metadata preserved in outputs

The chain spec is positional: `K101E+E138K` with `A+B` applies K101E to chain A
and E138K to chain B. If the chain spec has only one chain (e.g. `A`), that
chain is applied to every mutation in the combo.

## Outputs

### Local Pipeline
- `data/generated/<drug>/<mutation_label>/` - Minimized structures per mutation
- `results/metrics_summary.csv` - All metrics (binding proxy, contacts, H-bonds, pocket volume)
- `results/plots/<drug>_delta_metrics.png` - Bar charts of MUT-WT deltas

### Cluster FEP Pipeline
- `data/prepared/dor_4ncg/` - Minimized PDB structures
- `results/fep_runs/` - JSON results from cluster jobs
- `results/ddg_summary.csv` - ΔΔG values per mutation
- `results/correlation_analysis.csv` - Correlation statistics
- `results/plots/ddg_vs_fold_reduction.png` - Correlation plots

## CLI Reference

```bash
# Validation only (no OpenMM)
uv run python -m src.main --validate-only

# Verify mutations without minimization
uv run python -m src.main --verify-only

# Full local pipeline with replicates
uv run python -m src.main --replicates 3 --seed 42 --jitter-angstrom 0.1

# Force recomputation
uv run python -m src.main --force

# Cluster workflow
uv run python -m src.main --prepare-local --replicates 3 --seed 42 --jitter-angstrom 0.1
uv run python -m src.main --generate-slurm
uv run python -m src.main --collect-results

# SLURM customization
uv run python -m src.main --generate-slurm \
    --slurm-partition gpu \
    --slurm-time 4:00:00 \
    --slurm-memory 16G

# FEP parameters
uv run python -m src.main --prepare-local \
    --alchemy-equil-steps 10000 \
    --alchemy-prod-steps 25000 \
    --alchemy-sample-interval 200
```

## Dependencies

- **Structure processing**: pdbfixer, gemmi
- **Simulation**: openmm, openmmtools, openmmforcefields
- **Force fields**: openff-interchange, openff-forcefields, openff-toolkit
- **Analysis**: mdanalysis, rdkit
- **Core**: numpy, pandas, matplotlib
- **Utilities**: lxml, xmltodict, networkx, cachetools

## Technical Details

### Minimization Procedure

1. Load CIF structure with ligand from SDF
2. Apply AMBER14 protein (`ff14SB`) + DNA (`bsc1`) force fields
3. Use SMIRNOFF template for ligand (Gasteiger charges)
4. Apply harmonic restraints (500 kJ/mol/nm²) to atoms >8Å from ligand
5. Run `minimizeEnergy()` with restraints
6. Run second unrestrained minimization
7. Output minimized PDB

### Alchemical FEP Protocol

- **Lambda schedule**: 13 windows (1.0 → 0.0)
- **Equilibration**: 10,000 steps per window
- **Production**: 25,000 steps per window
- **Sample interval**: 200 steps
- **Free energy estimator**: Bennett Acceptance Ratio (BAR)
- **Runtime**: ~1 hour per leg on V100/A100

### Binding ΔG Calculation

```
ΔG_binding = ΔG_complex - ΔG_solvent
ΔΔG = ΔG_binding(mutant) - ΔG_binding(WT)
```

Positive ΔΔG indicates reduced binding affinity in the mutant (resistance).

### Structural Metrics

- **Contacts**: Atom pairs within 4.0Å between ligand and protein
- **H-bonds**: MDAnalysis HydrogenBondAnalysis (3.5Å, 135° cutoff)
- **Pocket volume**: Grid-based void volume within 8Å of ligand centroid

## Notes

- DNA is retained and lightly restrained during minimization
- The cluster workflow uses CUDA for GPU acceleration
- Local workflow defaults to Metal (macOS) or CPU
- Ligand SDFs are auto-generated from CIF metadata with RDKit hydrogens
