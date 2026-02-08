from __future__ import annotations

import argparse
import logging
from pathlib import Path

import pandas as pd

from .cluster import generate_slurm_script, run_result_collection
from .config import dor_4ncg_spec
from .dor_alchemy_pipeline import prepare_local_openmm_only_for_cluster
from .openmm.alchemy import AlchemicalConfig
from .plotting import plot_ddg_vs_fold_reduction, plot_lambda_profiles
from .utils import ensure_dirs, project_paths


def main(argv: list[str] | None = None) -> None:
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s"
    )
    parser = argparse.ArgumentParser(
        description="NNRTI OpenMM-only Sherlock pipeline"
    )
    parser.add_argument(
        "--prepare-local-openmm-only",
        action="store_true",
        help=(
            "Prepare minimized structures + prebuilt alchemical assets locally so "
            "Sherlock can run with OpenMM-only runtime dependencies."
        ),
    )
    parser.add_argument(
        "--generate-slurm",
        action="store_true",
        help="Generate SLURM submission script for Sherlock cluster.",
    )
    parser.add_argument(
        "--collect-results",
        action="store_true",
        help="Aggregate FEP results from cluster and compute correlations.",
    )
    parser.add_argument(
        "--replicates",
        type=int,
        default=1,
        help="Number of independent replicates per mutation (default: 1).",
    )
    parser.add_argument(
        "--jitter-angstrom",
        type=float,
        default=0.1,
        help="Random coordinate jitter (angstrom) applied before minimization.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Base RNG seed for jitter (default: None).",
    )
    parser.add_argument(
        "--mutation",
        type=str,
        default=None,
        help=(
            "Optional mutation label filter for preparation (e.g., V106A). "
            "WT is always included."
        ),
    )
    parser.add_argument(
        "--susceptibility-xlsx",
        type=Path,
        default=None,
        help="Path to DOR susceptibility workbook (default: data/DRM-susceptibilities.csv.xlsx).",
    )
    parser.add_argument(
        "--prepared-dir",
        type=Path,
        default=None,
        help="Directory for prepared CIFs (default: data/prepared/dor_4ncg).",
    )
    parser.add_argument(
        "--fep-manifest",
        type=Path,
        default=None,
        help="Path for FEP manifest CSV (default: results/fep_manifest.csv).",
    )
    parser.add_argument(
        "--slurm-script",
        type=Path,
        default=None,
        help="Path for generated SLURM script (default: scripts/sherlock/submit_fep.sh).",
    )
    parser.add_argument(
        "--fep-results-dir",
        type=Path,
        default=None,
        help="Directory containing FEP result JSONs (default: results/fep_runs).",
    )
    parser.add_argument(
        "--alchemy-equil-steps",
        type=int,
        default=10_000,
        help="Equilibration steps per lambda window for alchemical runs.",
    )
    parser.add_argument(
        "--alchemy-prod-steps",
        type=int,
        default=25_000,
        help="Production steps per lambda window for alchemical runs.",
    )
    parser.add_argument(
        "--alchemy-sample-interval",
        type=int,
        default=200,
        help="Sample interval (steps) for neighboring-window energy evaluations.",
    )
    parser.add_argument(
        "--trajectory-interval",
        type=int,
        default=2000,
        help="Step interval for trajectory frame writing on cluster.",
    )
    parser.add_argument(
        "--no-save-trajectories",
        action="store_true",
        help="Disable DCD trajectory output in generated SLURM worker commands.",
    )
    parser.add_argument(
        "--slurm-partition",
        type=str,
        default="gpu",
        help="SLURM partition for jobs (default: gpu).",
    )
    parser.add_argument(
        "--slurm-time",
        type=str,
        default="4:00:00",
        help="SLURM job time limit (default: 4:00:00).",
    )
    parser.add_argument(
        "--slurm-memory",
        type=str,
        default="16G",
        help="SLURM memory allocation (default: 16G).",
    )
    parser.add_argument(
        "--conda-env",
        type=str,
        default=None,
        help="Conda environment name to activate on cluster (e.g., nnrti).",
    )
    parser.add_argument(
        "--use-openmm-module",
        action="store_true",
        help="Generate SLURM script to use Sherlock chemistry/py-openmm module stack.",
    )
    args = parser.parse_args(argv)

    root = Path(__file__).resolve().parents[1]
    paths = project_paths(root)
    ensure_dirs([paths.generated, paths.results, paths.plots])

    susceptibility_xlsx = (
        args.susceptibility_xlsx
        if args.susceptibility_xlsx is not None
        else paths.data / "DRM-susceptibilities.csv.xlsx"
    )
    prepared_dir = (
        args.prepared_dir
        if args.prepared_dir is not None
        else paths.data / "prepared" / "dor_4ncg"
    )
    fep_manifest = (
        args.fep_manifest
        if args.fep_manifest is not None
        else paths.results / "fep_manifest.csv"
    )
    slurm_script = (
        args.slurm_script
        if args.slurm_script is not None
        else root / "scripts" / "sherlock" / "submit_fep.sh"
    )
    fep_results_dir = (
        args.fep_results_dir
        if args.fep_results_dir is not None
        else paths.results / "fep_runs"
    )
    alchemy_cfg = AlchemicalConfig(
        equilibration_steps=args.alchemy_equil_steps,
        production_steps=args.alchemy_prod_steps,
        sample_interval=args.alchemy_sample_interval,
    )

    if args.prepare_local_openmm_only:
        selected_mutations = {args.mutation} if args.mutation else None
        fep_tasks = prepare_local_openmm_only_for_cluster(
            root=root,
            susceptibility_xlsx=susceptibility_xlsx,
            prepared_dir=prepared_dir,
            fep_manifest_path=fep_manifest,
            fep_results_dir=fep_results_dir,
            replicates=args.replicates,
            jitter_seed_base=args.seed,
            jitter_angstrom=args.jitter_angstrom,
            alchemy_config=alchemy_cfg,
            selected_mutations=selected_mutations,
        )
        n_structures = len({(task.mutation, task.replicate) for task in fep_tasks})
        logging.info(
            "Prepared %d structure/replicate setups and %d FEP tasks",
            n_structures,
            len(fep_tasks),
        )
        logging.info("FEP manifest: %s", fep_manifest)
        return

    if args.generate_slurm:
        if not fep_manifest.exists():
            logging.error("FEP manifest not found: %s", fep_manifest)
            logging.error(
                "Run --prepare-local-openmm-only first to generate the manifest."
            )
            return
        output = generate_slurm_script(
            manifest_path=fep_manifest,
            output_script=slurm_script,
            partition=args.slurm_partition,
            time_limit=args.slurm_time,
            memory=args.slurm_memory,
            equil_steps=args.alchemy_equil_steps,
            prod_steps=args.alchemy_prod_steps,
            sample_interval=args.alchemy_sample_interval,
            trajectory_interval=args.trajectory_interval,
            save_trajectories=not args.no_save_trajectories,
            conda_env=args.conda_env,
            use_openmm_module=args.use_openmm_module,
        )
        logging.info("Generated SLURM script: %s", output)
        return

    if args.collect_results:
        if not fep_manifest.exists():
            logging.error("FEP manifest not found: %s", fep_manifest)
            return
        spec = dor_4ncg_spec(root)
        _, _, ddg_df = run_result_collection(
            manifest_path=fep_manifest,
            fep_results_dir=fep_results_dir,
            output_dir=paths.results,
            ligand_resname=spec.structure.ligand_resname,
        )
        logging.info("Wrote %s", paths.results / "ddg_summary.csv")
        logging.info("Wrote %s", paths.results / "correlation_analysis.csv")
        try:
            plot_ddg_vs_fold_reduction(ddg_df, paths)
            logging.info("Wrote %s", paths.plots / "ddg_vs_fold_reduction.png")
        except Exception as exc:
            logging.warning("Could not generate plot: %s", exc)
        try:
            lambda_summary_path = paths.results / "lambda_window_summary.csv"
            if lambda_summary_path.exists():
                lambda_summary_df = pd.read_csv(lambda_summary_path)
                plot_lambda_profiles(lambda_summary_df, paths)
                logging.info("Wrote %s", paths.plots / "lambda_profile_complex.png")
                logging.info("Wrote %s", paths.plots / "lambda_profile_solvent.png")
        except Exception as exc:
            logging.warning("Could not generate lambda-profile plots: %s", exc)
        return

    parser.error(
        "No action requested. Use one of: --prepare-local-openmm-only, "
        "--generate-slurm, --collect-results"
    )


if __name__ == "__main__":
    main()
