# NNRTI Resistance Mechanisms

Understanding structural mechanisms of Non-Nucleoside Reverse Transcriptase Inhibitor (NNRTI) resistance in HIV-1.

This repository computes binding free energy changes (ΔΔG) for drug resistance mutations using alchemical free energy perturbation (FEP) on GPU clusters.

## Quick Start: Cluster FEP Pipeline

### Overview

```
┌─────────────────────────────────────────────────────────────────────────────┐
│ LOCAL (CPU)                                                                 │
│                                                                             │
│   CIF structure ──► Mutagenesis ──► Min/solvate ──► FEP Manifest (CSV)      │
│                     (PDBFixer)       (assets)         (N tasks)            │
└─────────────────────────────────────────────────────────────────────────────┘
                                        │
                                        ▼ rsync to cluster
┌─────────────────────────────────────────────────────────────────────────────┐
│ CLUSTER (GPU)                                                               │
│                                                                             │
│   For each task (mutation × replicate × leg):                               │
│     1. Load prebuilt alchemical system                                     │
│     2. Run FEP: annihilate ligand across 13 λ windows                      │
│     3. Save ΔG to JSON                                                     │
└─────────────────────────────────────────────────────────────────────────────┘
                                        │
                                        ▼ rsync results back
┌─────────────────────────────────────────────────────────────────────────────┐
│ LOCAL (CPU)                                                                 │
│                                                                             │
│   Collect JSONs ──► Compute ΔΔG ──► Structural metrics ──► Correlations    │
│                     (ΔG_mut - ΔG_wt)  (contacts, H-bonds)   (vs resistance) │
└─────────────────────────────────────────────────────────────────────────────┘
```

### Step 1: Prepare Locally (recommended: prebuild OpenMM-only assets)

```bash
# Pilot run: WT + V106A only (OpenMM-only runtime on Sherlock)
# Use a conda environment with OpenMM + OpenMMTools for preparation.
conda activate nnrti-prep
export OPENMM_PLATFORM=CPU
python -m src.main \
  --prepare-local-openmm-only \
  --mutation V106A \
  --replicates 1 \
  --seed 42 \
  --jitter-angstrom 0.1
```

**Outputs:**
- `data/prepared/dor_4ncg/wt_4ncg.cif` - Wild-type structure
- `data/prepared/dor_4ncg/mut_v106a.cif` - V106A mutant structure
- `results/fep_runs/*/rep_*/assets/*_system.xml` - prebuilt alchemical systems
- `results/fep_runs/*/rep_*/assets/*_start.pdb` - minimized solvated topologies
- `results/fep_manifest.csv` - Task manifest (WT + V106A, complex/solvent legs)

**What's in the manifest:**
| Column | Description |
|--------|-------------|
| `task_id` | 0..N, used by SLURM array |
| `mutation` | "WT" or mutation label (e.g., "V106A") |
| `replicate` | 1, 2, or 3 |
| `leg` | "complex" or "solvent" (see below) |
| `input_cif` | Source CIF path used during local preparation |
| `jitter_seed` | Random seed for coordinate perturbation |
| `prepared_system_xml` | Prebuilt alchemical OpenMM system |
| `prepared_topology_pdb` | Starting solvated topology/coordinates |

### Step 2: Generate SLURM Script

```bash
python -m src.main --generate-slurm --use-openmm-module
```

**Output:** `scripts/sherlock/submit_fep.sh`

### Step 3: Setup on Sherlock (one-time)

For this OpenMM-only Sherlock path, use the site OpenMM module.

```bash
# SSH to Sherlock (replace <sunet-id> with your SUNet ID)
ssh <sunet-id>@login.sherlock.stanford.edu

# Load Sherlock OpenMM stack
ml chemistry py-openmm/8.1.1_py312

# Quick check (use python3 on Sherlock)
python3 -c "import openmm; print(openmm.__version__)"
```

If module import fails, fall back to a minimal conda environment:

```bash
ml miniforge/24.11.0-0
mamba create -n nnrti python=3.12 -y
mamba activate nnrti

# Runtime deps needed on Sherlock for prebuilt-asset mode
conda install -c conda-forge openmm numpy -y
```

Note: `openmmtools`, `openff*`, `openmmforcefields`, and `pdbfixer` are needed locally for preparation, not on Sherlock in this mode.

### Step 4: Transfer to Cluster

```bash
# Replace <sunet-id> with your SUNet ID
rsync -avz --exclude='.venv' --exclude='.git' \
    . <sunet-id>@login.sherlock.stanford.edu:/scratch/users/<sunet-id>/nnrti-mechanisms/

# For OpenMM-only assets, sync the prepared assets directory as well
rsync -avz results/fep_runs/ <sunet-id>@login.sherlock.stanford.edu:/scratch/users/<sunet-id>/nnrti-mechanisms/results/fep_runs/
```

### Step 5: Submit Jobs on Cluster

```bash
# On Sherlock (after ssh-ing in)
cd /scratch/users/<sunet-id>/nnrti-mechanisms
sbatch scripts/sherlock/submit_fep.sh

# Monitor progress
squeue -u <sunet-id>
```

For WT + V106A with 1 replicate, this submits 4 GPU jobs (2 mutations × 1 replicate × 2 legs).

**GPU allocation note:** The default `sh_dev -g 1` on Sherlock may land you
on an NVIDIA A30 with **MIG (Multi-Instance GPU) enabled**, giving you only
a ~6 GB / 14 SM slice instead of the full 24 GB / 56 SM GPU. This makes the
445K-atom complex leg extremely slow. For interactive testing, use `salloc`
on the `gpu` partition directly:

```bash
salloc -p gpu --gres=gpu:1 --time=4:00:00 --mem=32G
```

Verify you have a full GPU (no MIG section) with `nvidia-smi` after
allocation.

### Step 5b: Fix manifest paths on Sherlock

Local manifests store absolute paths from your workstation. Rewrite them on Sherlock before submitting jobs:

```bash
python3 scripts/sherlock/rewrite_manifest_paths.py
```

### Step 6: Transfer Results Back

```bash
# Replace <sunet-id> with your SUNet ID
rsync -avz <sunet-id>@login.sherlock.stanford.edu:/scratch/users/<sunet-id>/nnrti-mechanisms/results/fep_runs/ \
    results/fep_runs/
```

### Step 7: Collect and Analyze

```bash
python -m src.main --collect-results
```

**Outputs:**
- `results/ddg_summary.csv` - ΔΔG per mutation (mean ± std across replicates)
- `results/structural_metrics.csv` - Ensemble-averaged contacts, H-bonds, pocket volume
- `results/lambda_window_profiles.csv` - Per-window ΔG profile for each task
- `results/lambda_window_summary.csv` - WT/mutant window profile summary across replicates
- `results/correlation_analysis.csv` - Correlations vs fold-reduction for ΔΔG and structural metrics
- `results/plots/ddg_vs_fold_reduction.png`
- `results/plots/lambda_profile_complex.png`
- `results/plots/lambda_profile_solvent.png`

---

## Understanding the FEP Calculation

### Why Two Legs?

Binding free energy is computed via a thermodynamic cycle:

```
ΔG_binding = ΔG_complex - ΔG_solvent
```

| Leg | System | What happens |
|-----|--------|--------------|
| **complex** | Protein + Ligand + Water | Ligand is "turned off" while bound to protein |
| **solvent** | Ligand + Water (no protein) | Ligand is "turned off" in bulk water |

The difference gives the free energy of transferring the ligand from water to the binding site.

### Lambda Schedule

Each leg runs 13 windows where the ligand is gradually decoupled:

```
λ = 1.0 → 0.95 → 0.9 → 0.8 → 0.7 → 0.6 → 0.5 → 0.4 → 0.3 → 0.2 → 0.1 → 0.05 → 0.0
    ▲                                                                              ▲
    │                                                                              │
  Ligand fully                                                            Ligand fully
  interacting                                                             decoupled
```

Free energy differences between adjacent windows are computed using BAR (Bennett Acceptance Ratio).

### ΔΔG Interpretation

```
ΔΔG = ΔG_binding(mutant) - ΔG_binding(WT)
```

| ΔΔG | Meaning |
|-----|---------|
| Positive | Mutation weakens binding → resistance |
| Negative | Mutation strengthens binding → sensitization |
| ~0 | No significant effect |

---

## CLI Reference

### Cluster Workflow Commands

```bash
# Step 1: Prepare OpenMM-only assets locally (example: WT + V106A, 1 replicate)
python -m src.main --prepare-local-openmm-only \
    --mutation V106A \
    --replicates 1 \
    --seed 42 \
    --jitter-angstrom 0.1

# Step 2: Generate SLURM script (OpenMM module on Sherlock)
python -m src.main --generate-slurm \
    --use-openmm-module \
    --slurm-partition gpu \
    --slurm-time 4:00:00 \
    --slurm-memory 16G

# Step 7: Collect results
python -m src.main --collect-results
```

### FEP Parameters

```bash
--alchemy-equil-steps 10000    # Equilibration steps per λ window
--alchemy-prod-steps 25000     # Production steps per λ window
--alchemy-sample-interval 200  # Steps between energy samples
--trajectory-interval 2000     # Steps between saved DCD frames
```

### All Options

| Flag | Default | Description |
|------|---------|-------------|
| `--prepare-local-openmm-only` | - | Prebuild OpenMM-only assets for cluster |
| `--mutation` | None | Filter to a specific mutation label (e.g., V106A) |
| `--generate-slurm` | - | Generate SLURM submission script |
| `--collect-results` | - | Aggregate FEP results and compute ΔΔG |
| `--replicates` | 1 | Number of independent replicates |
| `--seed` | None | Base random seed for jitter |
| `--jitter-angstrom` | 0.1 | Coordinate perturbation magnitude |
| `--trajectory-interval` | 2000 | Steps between saved trajectory frames |
| `--no-save-trajectories` | False | Disable DCD trajectory output during FEP |
| `--conda-env` | None | Conda environment name on cluster |
| `--use-openmm-module` | False | Generate SLURM script for Sherlock OpenMM module |
| `--slurm-partition` | gpu | SLURM partition |
| `--slurm-time` | 4:00:00 | Job time limit |
| `--slurm-memory` | 16G | Memory per task |

---

## Input Data

### Structures
- `data/structures/4NCG.cif` - RT/DOR crystal structure (used for FEP)

### Susceptibility Data
- `data/DRM-susceptibilities.csv.xlsx` - DOR fold-reduction values from literature

---

## Output Files

### After `--prepare-local-openmm-only`
```
data/prepared/dor_4ncg/
├── wt_4ncg.cif              # Wild-type structure
├── mut_V106A.cif            # Single mutations
├── mut_V106M.cif
├── mut_L100I_K103N.cif      # Double mutations
└── ...

results/
└── fep_manifest.csv         # task definitions (absolute paths)
```

### After Cluster Run
```
results/fep_runs/
├── wt/
│   ├── rep_01/
│   │   ├── wt_minimized_rep01.pdb
│   │   ├── wt_complex_rep01.json
│   │   ├── wt_complex_rep01.dcd
│   │   ├── wt_complex_rep01_physical_lambda1.dcd
│   │   ├── wt_solvent_rep01.json
│   │   ├── wt_solvent_rep01.dcd
│   │   └── assets/
│   │       ├── wt_complex_rep01_start.pdb
│   │       ├── wt_complex_rep01_system.xml
│   │       ├── wt_solvent_rep01_start.pdb
│   │       └── wt_solvent_rep01_system.xml
│   ├── rep_02/
│   └── rep_03/
├── V106A/
│   ├── rep_01/
│   │   ├── V106A_minimized_rep01.pdb
│   │   ├── V106A_complex_rep01.json
│   │   ├── V106A_complex_rep01.dcd
│   │   ├── V106A_complex_rep01_physical_lambda1.dcd
│   │   ├── V106A_solvent_rep01.json
│   │   ├── V106A_solvent_rep01.dcd
│   │   └── assets/
│   │       ├── V106A_complex_rep01_start.pdb
│   │       ├── V106A_complex_rep01_system.xml
│   │       ├── V106A_solvent_rep01_start.pdb
│   │       └── V106A_solvent_rep01_system.xml
│   └── ...
└── ...
```

### After `--collect-results`
```
results/
├── ddg_summary.csv          # ΔΔG per mutation
├── ddg_full.csv             # All replicates
├── structural_metrics.csv   # Ensemble contacts, H-bonds, pocket volume
├── lambda_window_profiles.csv
├── lambda_window_summary.csv
├── correlation_analysis.csv # Correlations vs susceptibility
└── plots/
    ├── ddg_vs_fold_reduction.png
    ├── lambda_profile_complex.png
    └── lambda_profile_solvent.png
```

---

## Dependencies

```
Cluster runtime (FEP worker):
openmm numpy

Local prep / analysis tooling:
openmmtools openmmforcefields
openff-toolkit openff-forcefields
pdbfixer pandas openpyxl
gemmi rdkit mdanalysis matplotlib
```

---

## Technical Details

### Minimization (local prep)
- Force field: AMBER14 protein + DNA, SMIRNOFF ligand
- Restraints: 500 kJ/mol/nm² on atoms >8Å from ligand
- Two-stage: restrained then unrestrained
- No explicit solvent (gas phase, NoCutoff)

### Solvation (local prep)
- Water model: TIP3P
- Box padding: 1.0 nm
- Ionic strength: 0.15 M (Na⁺/Cl⁻)

### FEP Protocol
- Lambda protocol: two-phase decoupling (electrostatics first, then sterics)
  - Phase 1: λ_elec 1.0→0.0 with λ_sterics=1.0 (4 windows)
  - Phase 2: λ_sterics 1.0→0.0 with λ_elec=0.0 (9 windows)
- Electrostatics: PME with 1.0 nm cutoff
- Integrator: Langevin (300 K, 1/ps friction, 2 fs timestep)
- Barostat: Monte Carlo (1 bar)
- Free energy estimator: BAR
- Per-window stabilization: 100-iter minimization + 500-step warmup at 0.5 fs

### Sherlock Benchmarks (RTX 3090, DOR/4NCG system)

Measured from pilot run (WT + V106A, 1 replicate, job 15246033):

| Leg | Atoms | Wall time | Memory | GPU |
|-----|-------|-----------|--------|-----|
| Complex | 445,756 | ~1h 55m | 3.3 GB | RTX 3090 |
| Solvent | 1,065 | ~4 min | 140 MB | RTX 3090 |

Resource recommendations for SLURM:
- **Memory:** 8 GB is sufficient (peak 3.3 GB for complex leg)
- **Time:** 2.5 hrs per mutation (both legs sequentially on same GPU)
- **GPU:** 1× full GPU (no MIG). Use `salloc -p gpu` not `sh_dev`
- **Strategy:** Run both legs (complex + solvent) in the same job to halve
  the number of GPU allocations. The solvent leg adds only ~4 min overhead.

### Pilot Results (V106A, 1 replicate)

| Mutation | Complex ΔG | Solvent ΔG | Binding ΔG | ΔΔG |
|----------|-----------|-----------|-----------|-----|
| WT | -48.78 kJ/mol | -57.44 kJ/mol | +8.66 kJ/mol | 0 (ref) |
| V106A | -50.36 kJ/mol | -56.01 kJ/mol | +5.65 kJ/mol | -3.01 kJ/mol |

Note: Single replicate with short sampling (10K equil / 25K prod per window).
More replicates needed for reliable ΔΔG estimates.
