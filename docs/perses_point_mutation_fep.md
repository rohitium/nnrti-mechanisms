# Perses Point-Mutation FEP Protocol

This workflow computes the Jorgensen-style mutation thermodynamic cycle for
DOR resistance:

```text
holo: WT RT-DOR -> mutant RT-DOR   ΔG_holo
apo:  WT RT     -> mutant RT       ΔG_apo

ΔΔG_bind(mut - WT) = ΔG_holo - ΔG_apo
```

Positive `ΔΔG_bind` indicates that the mutation makes DOR binding less favorable
relative to WT.

## Environment

Create a clean environment rather than modifying `nnrti-prep`:

```bash
conda env create -f envs/nnrti-fep.yml
conda activate nnrti-fep
```

Perses uses OpenEye for SDF handling in the protein-mutation setup path, so the
environment includes `openeye-toolkits`. Confirm licensing works before launch:

```bash
python -c "from openeye import oechem; print(oechem.OEChemIsLicensed())"
python -c "import perses, openmmtools; print(perses.__version__)"
```

## Pilot

Run a short pilot for one mutation/replicate on Sherlock:

```bash
MUTATION_ALLOWLIST=V106A REPLICATES=1 FEP_N_CYCLES=500 \
  ./scripts/sherlock/submit_perses_point_mutation_fep.sh
```

Inspect `logs/pfep_V106A_1.*` and the corresponding
`results/analysis/perses_point_mutation_fep/V106A/**/summary.json`.

## Production

After the pilot passes, run a single-mutation panel:

```bash
MUTATION_ALLOWLIST=V106A,V106I,Y181C,Y188L,G190A,G190E,Y318F \
REPLICATES=1,2,3 FEP_N_CYCLES=5000 \
  ./scripts/sherlock/submit_perses_point_mutation_fep.sh
```

Summarize against the current MM/GBSA table:

```bash
PYTHONPATH=. python -m nnrti.cli.summarize_perses_point_mutation_fep
```

## Scope

This runner is for single point mutations such as `V106A` or `Y188L`. Double
mutants require either sequential alchemical cycles or a separate multi-mutation
setup and should not be mixed into this single-point panel without an explicit
design decision.
