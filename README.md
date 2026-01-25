# nnrti-mechanisms

Structural mechanisms of NNRTI resistance in HIV-1 patients.

This repository contains code, data, and analysis for evaluating the structural
impact of clinically relevant NNRTI resistance mutations on Rilpivirine (RPV)
and Doravirine (DOR) using cryo-EM RT/DNA/drug structures.

## What this repo does
- Curates drug-resistance mutation (DRM) lists in `data/DRMs.csv` (with chain IDs)
- Generates ligand SDFs directly from the CIF structures (with explicit H)
- Performs in silico mutagenesis on cryo-EM structures
- Minimizes WT and mutant complexes with OpenMM
- Computes relative binding proxies and structural metrics:
  - Binding energy proxy (OpenMM potential energy delta)
  - Ligand-protein contacts
  - H-bond count
  - Binding pocket volume proxy (A^3)
- Produces per-drug bar charts of mutant minus WT deltas

## Key inputs
- Structures:
  - `data/structures/7Z2D.cif` (RT/DNA/RPV)
  - `data/structures/7Z2G.cif` (RT/DNA/DOR)
- DRM list:
  - `data/DRMs.csv` (includes a `chain` column; chain-to-subunit mapping is parsed from CIFs)

## DRM table format
`data/DRMs.csv` must include:
- `drug`: `RPV` or `DOR`
- `mutation`: one or more mutations, e.g. `K101E` or `K101E+E138K`
- `chain`: chain id(s) for each mutation, e.g. `A` or `A+B`
- `category`, `notes`: optional metadata preserved in outputs

The chain spec is positional: `K101E+E138K` with `A+B` applies K101E to chain A
and E138K to chain B. If the chain spec has only one chain (e.g. `A`), that
chain is applied to every mutation in the combo.

## Outputs
- Minimized structures and artifacts:
  - `data/generated/rpv/<mutation_label>/` (one folder per mutation)
  - `data/generated/dor/<mutation_label>/` (one folder per mutation)
  - `data/generated/<drug>/wt/` (WT minimized structures)
- Metrics table:
  - `results/metrics_summary.csv`
  - `results/metrics_summary.xlsx` (two sheets: `RPV`, `DOR`; WT rows stored once per drug; no `structure`/chain/subunit columns)
- Plots:
  - `results/plots/rpv_delta_metrics.png`
  - `results/plots/dor_delta_metrics.png`

## Dependencies
- openmm
- pdbfixer
- openmmforcefields
- openff-toolkit
- openff-units
- openff-utilities
- openff-interchange
- openff-forcefields
- mdanalysis
- rdkit
- gemmi (only required for `src.ligand_from_cif`)
- numpy, pandas, matplotlib
- lxml, xmltodict, networkx, cachetools, python-constraint

## How to run
1) Generate ligand SDFs from CIF metadata:
```bash
uv run python -m src.ligand_from_cif
```

2) Run the DRM batch pipeline (uses `data/DRMs.csv` for all mutations/combos):
```bash
uv run python -m src.main
```

The pipeline writes `results/metrics_summary.csv` and updates the per-drug plots
in `results/plots/`.

## Pipeline logic (step-by-step)
1) Load `data/DRMs.csv` and filter rows by drug (`RPV`, `DOR`).
2) Parse chain-to-subunit mapping from the CIF files via `_entity_name_com` /
   `_entity.pdbx_description` + `_struct_asym` to label chains as `p66` or `p51`.
3) For each drug:
   - Minimize the WT complex once (`wt_minimized.cif/.pdb`) and compute metrics.
   - For each DRM row, apply mutations (single or combo) to the specified chain(s),
     minimize, and compute metrics.
4) Write `results/metrics_summary.csv` with both WT and MUT rows for every metric.
5) Plot per-drug MUT–WT deltas for each metric to `results/plots/*.png`.

## Data preprocessing
- Ligand SDFs are generated from CIF metadata using `src/ligand_from_cif.py`,
  with explicit hydrogens added by RDKit.
- For OpenMM compatibility, the ligand in original cryo-EM CIF is replaced by the SDF ligand
  (same residue name) before minimization.
- Nonstandard DNA residue `OMC` is converted to `DC` (deoxycytidine) for force
  field compatibility.
- All existing hydrogens in the CIF are removed before re-adding hydrogens to
  avoid duplicates and to ensure consistent protonation/geometry under the
  force field (CIF hydrogens can be incomplete or inconsistently named).

Note: `OMC` (5-methylcytosine) is converted to `DC` to keep the DNA compatible
with the standard AMBER DNA force field. This removes the methyl group and is
an approximation; in these structures it is part of the DNA aptamer and is
typically distant from the NNRTI pocket, so the impact on ligand-proximal
metrics is expected to be small.

## Minimization procedure
Minimization is performed with OpenMM on each WT and mutant complex:
- Input: the processed cryo-EM CIF with the ligand added.
- Force field: AMBER14 protein (`ff14SB`) + AMBER14 DNA (`bsc1`), with a SMIRNOFF
  template for the ligand (Gasteiger charges).
- Nonbonded interactions use `NoCutoff`; bonds to hydrogen are constrained.
- A harmonic positional restraint is applied to heavy atoms farther than
  8 Å from any ligand atom (ligand-adjacent region is left flexible).
  Restraint strength is 500 kJ/mol/nm^2.
- Minimization uses a Langevin integrator (300 K, 1/ps friction, 2 fs step)
  and runs `Simulation.minimizeEnergy()` (no dynamics).
- Outputs: minimized `*.cif` and a corresponding `*.pdb` for metric calculations.

## Binding proxy definition and interpretation
For each minimized structure, the pipeline computes:
- `E_complex`: potential energy of the full system (protein + DNA + ligand)
- `E_receptor`: potential energy after removing the ligand (protein + DNA only)
- `E_ligand`: potential energy of the ligand alone (same coordinates; no re-minimization)

The binding proxy is:
```
E_binding_proxy = E_complex − E_receptor − E_ligand
```

This is an interaction-energy-like proxy, not a binding free energy. It reuses
one minimized configuration and does not relax the receptor or ligand after
removal.

The plotted delta is:
```
Δ = (E_binding_proxy)_MUT − (E_binding_proxy)_WT
```
Positive Δ indicates a less favorable binding proxy in the mutant; negative Δ
indicates a more favorable binding proxy.

## Contacts definition and interpretation
Contacts are computed on the minimized structure using MDAnalysis. The pipeline:
- Selects the ligand by residue name (`resname`).
- Selects the protein with `protein and not resname <ligand>` (DNA excluded).
- Computes all pairwise distances between ligand atoms and protein atoms.
- Counts the number of atom pairs within 4.0 Å.

This yields a raw contact count:
```
contact_count = number of (ligand atom, protein atom) pairs with distance < 4.0 Å
```

The plotted delta is:
```
Δ Contacts = contact_count_MUT − contact_count_WT
```

Interpretation:
- Positive Δ means the mutant has more close atom–atom contacts (tighter packing).
- Negative Δ means fewer close contacts (looser packing).
- It is a geometric proxy, not weighted by atom type or interaction strength.
- Hydrogens are included if present in the minimized structure, so absolute
  values can be sensitive to hydrogen placement.

## H-bonds definition and interpretation
H-bonds are computed on the minimized structure using MDAnalysis'
`HydrogenBondAnalysis`:
- Donors: protein or ligand.
- Acceptors: protein or ligand.
- Hydrogens: atoms named `H*`.
- Only protein–ligand pairs are counted.
- Distance cutoff: 3.5 Å.
- Angle cutoff: 135°.

This yields:
```
hbond_count = number of protein–ligand H-bonds meeting the criteria
```

The plotted delta is:
```
Δ H-bonds = hbond_count_MUT − hbond_count_WT
```

Interpretation:
- Positive Δ means more protein–ligand H-bonds in the mutant.
- Negative Δ means fewer protein–ligand H-bonds in the mutant.
- It depends on hydrogen placement and the geometric criteria above.

## Binding pocket volume proxy definition and interpretation
The pocket volume proxy is computed on a cubic grid centered at the ligand:
- All grid points within an 8 Å radius of the ligand centroid are considered.
- Grid spacing is 0.5 Å (voxel volume = 0.125 Å^3).
- Receptor atoms are all non-hydrogen atoms excluding the ligand.
- A grid point is considered "free" if it is outside the van der Waals radius
  (element-specific) of every receptor atom.

The proxy is:
```
pocket_volume_proxy = number of free grid points * voxel volume (Å^3)
```

The plotted delta is:
```
Δ Pocket Volume = pocket_volume_proxy_MUT − pocket_volume_proxy_WT
```

Interpretation:
- Positive Δ means a larger empty pocket around the ligand in the mutant.
- Negative Δ means a more constricted pocket.
- This is a coarse geometric proxy, not a solvent-accessible or physical volume.

## Metrics table columns
Each row in `results/metrics_summary.csv` includes:
- `structure`: `RPV` or `DOR`
- `mutation`: original mutation string from the DRM table
- `chain`: chain id(s) used (e.g. `A` or `A+B`)
- `subunit`: derived subunit label(s) (e.g. `p66` or `p66+p51`)
- `category`, `notes`: passthrough from DRM table
- `state`: `WT` or `MUT`
- `metric`: `binding_proxy_kj_mol`, `contact_count`, `hbond_count`, `pocket_volume_proxy`
- `value`: metric value (binding energy in kJ/mol; pocket volume in A^3)

## Performance
- Mutations are processed in parallel across all CPU cores using multiprocessing.
- The WT minimization is performed once per drug and reused for all deltas.

## Notes
- DNA is retained in the complex and lightly restrained outside a ligand-centric
  shell during minimization to reduce drift.
- Energies are OpenMM potential energies used as relative binding proxies.
- The OpenMM build in the current venv uses CPU; Metal GPU support requires a
  conda-forge OpenMM build (not yet configured here).
