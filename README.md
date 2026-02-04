# NNRTI Resistance Mechanisms

Understanding structural mechanisms of Non-Nucleoside Reverse Transcriptase Inhibitor (NNRTI) resistance in HIV-1.

This repository computes binding free energy changes (ΔΔG) for drug resistance mutations using alchemical free energy perturbation (FEP) on GPU clusters.

## Quick Start: Cluster FEP Pipeline

### Overview

```
┌─────────────────────────────────────────────────────────────────────────────┐
│ LOCAL (CPU)                                                                 │
│                                                                             │
│   CIF structure ──► Mutagenesis ──► Mutant CIFs ──► FEP Manifest (CSV)     │
│                     (PDBFixer)       (14 structures)  (84 tasks)            │
└─────────────────────────────────────────────────────────────────────────────┘
                                        │
                                        ▼ rsync to cluster
┌─────────────────────────────────────────────────────────────────────────────┐
│ CLUSTER (GPU)                                                               │
│                                                                             │
│   For each task (mutation × replicate × leg):                               │
│     1. Minimize structure (OpenMM CUDA)                                     │
│     2. Solvate with TIP3P water + 0.15M ions                               │
│     3. Run FEP: annihilate ligand across 13 λ windows                      │
│     4. Save ΔG to JSON                                                      │
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

### Step 1: Prepare Locally (creates mutant CIFs + manifest)

```bash
uv run python -m src.main --prepare-local --replicates 3 --seed 42 --jitter-angstrom 0.1
```

**Outputs:**
- `data/prepared/dor_4ncg/wt_4ncg.cif` - Wild-type structure
- `data/prepared/dor_4ncg/mut_*.cif` - Mutant structures (13 mutations)
- `results/fep_manifest.csv` - Task manifest (84 tasks)

**What's in the manifest:**
| Column | Description |
|--------|-------------|
| `task_id` | 0-83, used by SLURM array |
| `mutation` | "WT" or mutation label (e.g., "V106A") |
| `replicate` | 1, 2, or 3 |
| `leg` | "complex" or "solvent" (see below) |
| `input_cif` | Path to CIF for minimization |
| `jitter_seed` | Random seed for coordinate perturbation |

### Step 2: Generate SLURM Script

```bash
uv run python -m src.main --generate-slurm --conda-env nnrti
```

**Output:** `scripts/sherlock/submit_fep.sh`

### Step 3: Setup on Sherlock (one-time)

Some packages (openmmtools, pdbfixer) are only available via conda-forge.

```bash
# SSH to Sherlock (replace <sunet-id> with your SUNet ID)
ssh <sunet-id>@login.sherlock.stanford.edu

# Check available conda modules
module spider conda
module avail anaconda

# Load conda module (name may vary - check output above)
ml anaconda3  # or whatever is available

# Create conda environment (one-time)
conda create -n nnrti python=3.12 -y
conda activate nnrti

# Install dependencies from conda-forge
conda install -c conda-forge \
    openmm openmmtools openmmforcefields \
    openff-toolkit openff-forcefields \
    pdbfixer gemmi rdkit mdanalysis \
    numpy pandas matplotlib -y
```

Note: Installation may take a while due to dependency resolution.

### Step 4: Transfer to Cluster

```bash
# Replace <sunet-id> with your SUNet ID
rsync -avz --exclude='.venv' --exclude='.git' \
    . <sunet-id>@login.sherlock.stanford.edu:/scratch/users/<sunet-id>/nnrti-mechanisms/
```

### Step 5: Submit Jobs on Cluster

```bash
# On Sherlock (after ssh-ing in)
cd /scratch/users/<sunet-id>/nnrti-mechanisms
sbatch scripts/sherlock/submit_fep.sh

# Monitor progress
squeue -u <sunet-id>
```

This submits 84 parallel GPU jobs (one per task).

### Step 6: Transfer Results Back

```bash
# Replace <sunet-id> with your SUNet ID
rsync -avz <sunet-id>@login.sherlock.stanford.edu:/scratch/users/<sunet-id>/nnrti-mechanisms/results/fep_runs/ \
    results/fep_runs/
```

### Step 7: Collect and Analyze

```bash
uv run python -m src.main --collect-results
```

**Outputs:**
- `results/ddg_summary.csv` - ΔΔG per mutation (mean ± std across replicates)
- `results/structural_metrics.csv` - Contacts, H-bonds, pocket volume
- `results/correlation_analysis.csv` - Pearson/Spearman vs fold-reduction
- `results/plots/ddg_vs_fold_reduction.png`

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
# Step 1: Prepare mutant structures and manifest
uv run python -m src.main --prepare-local \
    --replicates 3 \
    --seed 42 \
    --jitter-angstrom 0.1

# Step 2: Generate SLURM script
uv run python -m src.main --generate-slurm \
    --conda-env nnrti \
    --slurm-partition gpu \
    --slurm-time 4:00:00 \
    --slurm-memory 16G

# Step 7: Collect results
uv run python -m src.main --collect-results
```

### FEP Parameters

```bash
--alchemy-equil-steps 10000    # Equilibration steps per λ window
--alchemy-prod-steps 25000     # Production steps per λ window
--alchemy-sample-interval 200  # Steps between energy samples
```

### All Options

| Flag | Default | Description |
|------|---------|-------------|
| `--prepare-local` | - | Create mutant CIFs and FEP manifest |
| `--generate-slurm` | - | Generate SLURM submission script |
| `--collect-results` | - | Aggregate FEP results and compute ΔΔG |
| `--replicates` | 1 | Number of independent replicates |
| `--seed` | None | Base random seed for jitter |
| `--jitter-angstrom` | 0.0 | Coordinate perturbation magnitude |
| `--conda-env` | None | Conda environment name on cluster |
| `--slurm-partition` | gpu | SLURM partition |
| `--slurm-time` | 4:00:00 | Job time limit |
| `--slurm-memory` | 16G | Memory per task |

---

## Input Data

### Structures
- `data/structures/4NCG.cif` - RT/DOR crystal structure (used for FEP)
- `data/structures/7Z2D.cif` - RT/DNA/RPV cryo-EM structure
- `data/structures/7Z2G.cif` - RT/DNA/DOR cryo-EM structure

### Susceptibility Data
- `data/DRM-susceptibilities.csv.xlsx` - DOR fold-reduction values from literature

---

## Output Files

### After `--prepare-local`
```
data/prepared/dor_4ncg/
├── wt_4ncg.cif              # Wild-type structure
├── mut_V106A.cif            # Single mutations
├── mut_V106M.cif
├── mut_L100I_K103N.cif      # Double mutations
└── ...

results/
└── fep_manifest.csv         # 84 task definitions
```

### After Cluster Run
```
results/fep_runs/
├── wt/
│   ├── rep_01/
│   │   ├── wt_minimized_rep01.pdb
│   │   ├── wt_complex_rep01.json
│   │   └── wt_solvent_rep01.json
│   ├── rep_02/
│   └── rep_03/
├── V106A/
│   ├── rep_01/
│   │   ├── V106A_minimized_rep01.pdb
│   │   ├── V106A_complex_rep01.json
│   │   └── V106A_solvent_rep01.json
│   └── ...
└── ...
```

### After `--collect-results`
```
results/
├── ddg_summary.csv          # ΔΔG per mutation
├── ddg_full.csv             # All replicates
├── structural_metrics.csv   # Contacts, H-bonds, pocket volume
├── correlation_analysis.csv # Pearson/Spearman statistics
└── plots/
    └── ddg_vs_fold_reduction.png
```

---

## Dependencies

```
openmm openmmtools openmmforcefields
openff-toolkit openff-forcefields
pdbfixer gemmi rdkit mdanalysis
numpy pandas matplotlib
```

---

## Technical Details

### Minimization (on cluster)
- Force field: AMBER14 protein + DNA, SMIRNOFF ligand
- Restraints: 500 kJ/mol/nm² on atoms >8Å from ligand
- Two-stage: restrained then unrestrained
- No explicit solvent (gas phase, NoCutoff)

### Solvation (on cluster)
- Water model: TIP3P
- Box padding: 1.0 nm
- Ionic strength: 0.15 M (Na⁺/Cl⁻)

### FEP Protocol
- Electrostatics: PME with 1.0 nm cutoff
- Integrator: Langevin (300 K, 1/ps friction, 2 fs timestep)
- Barostat: Monte Carlo (1 bar)
- Free energy estimator: BAR

### Expected Runtime
- ~1 hour per leg on V100/A100
- Total: ~84 GPU-hours for full dataset
