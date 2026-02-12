"""
Snakemake workflow for NNRTI resistance mechanism simulations.

Runs on Sherlock login node; GPU MD rules are submitted via SLURM executor.
Usage:
    snakemake -n                                       # dry run
    snakemake --profile workflow/profiles/sherlock      # production run on Sherlock
    snakemake -j1                                      # local run (single core, for testing)
"""

from pathlib import Path

import pandas as pd

from src.analysis.susceptibility import load_dor_susceptibilities
from src.utils.mutations import deterministic_seed, sanitize_label

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
configfile: "workflow/config.yaml"

REPLICATES = list(range(1, config["replicates"] + 1))
REPS_STR = [f"{r:02d}" for r in REPLICATES]

# Build mutation list from susceptibility data at parse time.
_dor_df = load_dor_susceptibilities(
    Path(config["susceptibility_xlsx"]), default_chain="A"
)
_MUTANT_NAMES = _dor_df["mutation"].tolist()
MUTATIONS = ["WT"] + _MUTANT_NAMES

# Mappings between mutation name <-> safe label.
SAFE_LABELS = {m: (sanitize_label(m) if m != "WT" else "wt") for m in MUTATIONS}
LABEL_TO_MUTATION = {v: k for k, v in SAFE_LABELS.items()}
ALL_SAFE_LABELS = list(SAFE_LABELS.values())

# Fold-reduction lookup.
_FOLD_REDUCTION = dict(
    zip(_dor_df["mutation"].tolist(), _dor_df["dor_fold_reduction"].tolist())
)


def _cif_for_label(safe_label):
    """Return the CIF path for a given safe_label."""
    if safe_label == "wt":
        return "data/prepared/dor_4ncg/wt_4ncg.cif"
    return f"data/prepared/dor_4ncg/mut_{safe_label}.cif"


# ---------------------------------------------------------------------------
# Wildcard constraints
# ---------------------------------------------------------------------------
wildcard_constraints:
    safe_label="[A-Za-z0-9_]+",
    rep="[0-9]{2}",


# ---------------------------------------------------------------------------
# Target rule
# ---------------------------------------------------------------------------
rule all:
    input:
        "results/ddg_full.csv",
        "results/plots/all_metrics_vs_fold_reduction.png",
        "results/plots/boundness_qc_min_distance.png",
        "results/plots/rmsd_convergence.png",
        "results/plots/com_distance_convergence.png",


# ---------------------------------------------------------------------------
# Rule 1: Copy WT CIF
# ---------------------------------------------------------------------------
rule prep_wt_cif:
    input:
        config["wt_cif_source"],
    output:
        "data/prepared/dor_4ncg/wt_4ncg.cif",
    shell:
        "mkdir -p $(dirname {output}) && cp {input} {output}"


# ---------------------------------------------------------------------------
# Rule 2: Apply mutations to produce mutant CIF
# ---------------------------------------------------------------------------
rule prep_mutant_cif:
    input:
        wt_cif="data/prepared/dor_4ncg/wt_4ncg.cif",
        xlsx=config["susceptibility_xlsx"],
    output:
        "data/prepared/dor_4ncg/mut_{safe_label}.cif",
    params:
        mutation=lambda wc: LABEL_TO_MUTATION[wc.safe_label],
    resources:
        mem_mb=4000,
        runtime=10,
    script:
        "workflow/scripts/prep_mutant_cif.py"


# ---------------------------------------------------------------------------
# Rule 3: Minimize + solvate for one replicate
# ---------------------------------------------------------------------------
rule prep_replicate:
    input:
        cif=lambda wc: _cif_for_label(wc.safe_label),
        ligand_sdf=config["ligand_sdf"],
    output:
        minimized_pdb="results/md_runs/{safe_label}/rep_{rep}/{safe_label}_minimized_rep{rep}.pdb",
        system_xml="results/md_runs/{safe_label}/rep_{rep}/assets/{safe_label}_md_rep{rep}_system.xml",
        topology_pdb="results/md_runs/{safe_label}/rep_{rep}/assets/{safe_label}_md_rep{rep}_start.pdb",
    params:
        seed=lambda wc: deterministic_seed(config["seed"], wc.safe_label, int(wc.rep)),
        jitter=config["jitter_angstrom"],
        ligand_resname=config["ligand_resname"],
    envmodules:
        "chemistry",
        "py-openmm/8.1.1_py312",
    resources:
        mem_mb=16000,
        runtime=30,
    script:
        "workflow/scripts/prep_replicate.py"


# ---------------------------------------------------------------------------
# Rule 4: Heating + production MD (GPU)
# ---------------------------------------------------------------------------
rule run_md:
    input:
        system_xml="results/md_runs/{safe_label}/rep_{rep}/assets/{safe_label}_md_rep{rep}_system.xml",
        topology_pdb="results/md_runs/{safe_label}/rep_{rep}/assets/{safe_label}_md_rep{rep}_start.pdb",
    output:
        result_json="results/md_runs/{safe_label}/rep_{rep}/{safe_label}_rep{rep}.json",
        final_pdb="results/md_runs/{safe_label}/rep_{rep}/{safe_label}_rep{rep}_md_final.pdb",
        analysis_dcd="results/md_runs/{safe_label}/rep_{rep}/{safe_label}_rep{rep}_analysis.dcd",
        analysis_topo="results/md_runs/{safe_label}/rep_{rep}/{safe_label}_rep{rep}_analysis_topology.pdb",
        state_csv="results/md_runs/{safe_label}/rep_{rep}/{safe_label}_rep{rep}_md_state.csv",
    params:
        mutation=lambda wc: LABEL_TO_MUTATION[wc.safe_label],
        heating_ps=config["md"]["heating_ps"],
        production_ns=config["md"]["production_ns"],
        report_interval=config["md"]["report_interval"],
        checkpoint_interval=config["md"]["checkpoint_interval"],
        ligand_sdf=config["ligand_sdf"],
        ligand_resname=config["ligand_resname"],
        fold_reduction=lambda wc: _FOLD_REDUCTION.get(
            LABEL_TO_MUTATION[wc.safe_label]
        ),
    envmodules:
        "chemistry",
        "py-openmm/8.1.1_py312",
    resources:
        slurm_partition=config["slurm"]["partition"],
        slurm_extra="'--gres={}'".format(config["slurm"]["gres"]),
        mem_mb=config["slurm"]["mem_mb"],
        runtime=config["slurm"]["time_minutes"],
    script:
        "workflow/scripts/run_md.py"


# ---------------------------------------------------------------------------
# Rule 5: Collect results + MM/GBSA + structural metrics + aggregation
# ---------------------------------------------------------------------------
rule collect_and_analyze:
    input:
        jsons=expand(
            "results/md_runs/{safe_label}/rep_{rep}/{safe_label}_rep{rep}.json",
            zip,
            safe_label=[sl for sl in ALL_SAFE_LABELS for _ in REPS_STR],
            rep=[r for _ in ALL_SAFE_LABELS for r in REPS_STR],
        ),
        xlsx=config["susceptibility_xlsx"],
        ligand_sdf=config["ligand_sdf"],
    output:
        ddg_full="results/ddg_full.csv",
        correlation="results/correlation_analysis.csv",
        mmgbsa="results/mmgbsa_replicate_metrics.csv",
        structural="results/structural_metrics.csv",
        boundness="results/boundness_qc.csv",
        rmsd_profiles="results/rmsd_ca_profiles.csv",
        com_profiles="results/com_distance_profiles.csv",
    params:
        mmgbsa_snapshots=config["analysis"]["mmgbsa_snapshots"],
        mmgbsa_discard_fraction=config["analysis"]["mmgbsa_discard_fraction"],
        metric_frame_stride=config["analysis"]["metric_frame_stride"],
        metric_max_frames=config["analysis"]["metric_max_frames"],
        ligand_resname=config["ligand_resname"],
    envmodules:
        "chemistry",
        "py-openmm/8.1.1_py312",
    resources:
        mem_mb=32000,
        runtime=120,
    script:
        "workflow/scripts/collect_and_analyze.py"


# ---------------------------------------------------------------------------
# Rule 6: Generate all plots
# ---------------------------------------------------------------------------
rule generate_plots:
    input:
        ddg_full="results/ddg_full.csv",
        boundness="results/boundness_qc.csv",
        rmsd="results/rmsd_ca_profiles.csv",
        com="results/com_distance_profiles.csv",
    output:
        "results/plots/all_metrics_vs_fold_reduction.png",
        "results/plots/boundness_qc_min_distance.png",
        "results/plots/rmsd_convergence.png",
        "results/plots/com_distance_convergence.png",
    resources:
        mem_mb=4000,
        runtime=10,
    script:
        "workflow/scripts/generate_plots.py"
