from __future__ import annotations

import argparse
import logging
from pathlib import Path

import pandas as pd

from .cluster import generate_slurm_script, run_result_collection
from .config import dor_spec, rpv_spec
from .dor_alchemy_pipeline import (
    prepare_dor_local_structures,
    prepare_local_with_fep_manifest,
    run_dor_alchemical_manifest,
    summarize_dor_correlations,
)
from .drm_io import load_drms
from .metrics_io import write_metrics_xlsx
from .mutation.runner import run_mutations
from .numbering import detect_numbering_scheme
from .openmm.alchemy import AlchemicalConfig
from .plotting import plot_delta_metrics, plot_ddg_vs_fold_reduction
from .utils import ensure_dirs, load_chain_subunits, load_residue_mappings, project_paths
from .validation import validate_mutations, verify_mutations_only


def main(argv: list[str] | None = None) -> None:
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s"
    )
    parser = argparse.ArgumentParser(description="NNRTI DRM pipeline")
    parser.add_argument(
        "--validate-only",
        action="store_true",
        help="Validate DRM substitutions against the CIF sequences without OpenMM.",
    )
    parser.add_argument(
        "--verify-only",
        action="store_true",
        help="Apply DRMs and verify substitutions without OpenMM.",
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
        default=0.0,
        help="Random coordinate jitter (angstrom) applied before minimization.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Base RNG seed for jitter (default: None).",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Recompute metrics even if results/metrics_summary.csv exists.",
    )
    parser.add_argument(
        "--prepare-dor-local",
        action="store_true",
        help="Prepare DOR (4NCG) WT/mutant CIFs from susceptibility sheet for cluster runs.",
    )
    parser.add_argument(
        "--run-dor-manifest",
        action="store_true",
        help="Run DOR alchemical workflow from a prepared manifest.",
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
        "--manifest-csv",
        type=Path,
        default=None,
        help="Path for prepared manifest CSV (default: results/dor_4ncg_manifest.csv).",
    )
    parser.add_argument(
        "--mutation-index",
        type=int,
        default=None,
        help="Optional 0-based row index into manifest for single-mutation execution.",
    )
    parser.add_argument(
        "--alchemy-equil-steps",
        type=int,
        default=10000,
        help="Equilibration steps per lambda window for alchemical runs.",
    )
    parser.add_argument(
        "--alchemy-prod-steps",
        type=int,
        default=40000,
        help="Production steps per lambda window for alchemical runs.",
    )
    parser.add_argument(
        "--alchemy-sample-interval",
        type=int,
        default=200,
        help="Sample interval (steps) for neighboring-window energy evaluations.",
    )

    parser.add_argument(
        "--prepare-local",
        action="store_true",
        help="Prepare structures, compute structural metrics, and generate FEP manifest for cluster.",
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
    manifest_csv = (
        args.manifest_csv
        if args.manifest_csv is not None
        else paths.results / "dor_4ncg_manifest.csv"
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

    if args.prepare_local:
        logging.info("Preparing local structures and generating FEP manifest")
        structural_metrics, fep_tasks = prepare_local_with_fep_manifest(
            root=root,
            susceptibility_xlsx=susceptibility_xlsx,
            prepared_dir=prepared_dir,
            fep_manifest_path=fep_manifest,
            structural_metrics_path=paths.results / "structural_metrics.csv",
            fep_results_dir=fep_results_dir,
            replicates=args.replicates,
            jitter_seed_base=args.seed,
            jitter_angstrom=args.jitter_angstrom,
        )
        logging.info(
            "Prepared %d structures, generated %d FEP tasks",
            len(structural_metrics) // args.replicates,
            len(fep_tasks),
        )
        logging.info("Structural metrics: %s", paths.results / "structural_metrics.csv")
        logging.info("FEP manifest: %s", fep_manifest)
        return

    if args.generate_slurm:
        if not fep_manifest.exists():
            logging.error("FEP manifest not found: %s", fep_manifest)
            logging.error("Run --prepare-local first to generate the manifest.")
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
        )
        logging.info("Generated SLURM script: %s", output)
        return

    if args.collect_results:
        if not fep_manifest.exists():
            logging.error("FEP manifest not found: %s", fep_manifest)
            return
        structural_metrics_path = paths.results / "structural_metrics.csv"
        ddg_summary, correlations, ddg_df = run_result_collection(
            manifest_path=fep_manifest,
            fep_results_dir=fep_results_dir,
            structural_metrics_path=structural_metrics_path,
            output_dir=paths.results,
        )
        logging.info("ΔΔG summary: %s", paths.results / "ddg_summary.csv")
        logging.info("Correlations: %s", paths.results / "correlation_analysis.csv")

        try:
            plot_ddg_vs_fold_reduction(ddg_df, paths)
            logging.info("Plot: %s", paths.plots / "ddg_vs_fold_reduction.png")
        except Exception as e:
            logging.warning("Could not generate plot: %s", e)
        return

    if args.prepare_dor_local:
        manifest = prepare_dor_local_structures(
            root=root,
            susceptibility_xlsx=susceptibility_xlsx,
            prepared_dir=prepared_dir,
            manifest_csv=manifest_csv,
        )
        logging.info("Prepared %d DOR mutation structures at %s", len(manifest), prepared_dir)
        logging.info("Manifest written to %s", manifest_csv)
        return

    if args.run_dor_manifest:
        df = run_dor_alchemical_manifest(
            root=root,
            manifest_csv=manifest_csv,
            output_dir=paths.generated / "dor",
            replicates=args.replicates,
            jitter_seed_base=args.seed,
            jitter_angstrom=args.jitter_angstrom,
            alchemy_config=alchemy_cfg,
            mutation_index=args.mutation_index,
        )
        output_csv = paths.results / "metrics_summary.csv"
        df.to_csv(output_csv, index=False)
        write_metrics_xlsx(df, paths.results / "metrics_summary.xlsx")
        corr = summarize_dor_correlations(df)
        corr.to_csv(paths.results / "dor_correlations.csv", index=False)
        plot_delta_metrics(df, paths)
        logging.info("Wrote %s and %s", output_csv, paths.results / "dor_correlations.csv")
        return

    drms = load_drms(paths.data / "DRMs.csv")
    if args.validate_only:
        validate_mutations(drms, root)
        return
    if args.verify_only:
        verify_mutations_only(drms, root)
        return

    output_csv = paths.results / "metrics_summary.csv"
    if output_csv.exists() and not args.force:
        logging.info("Using existing metrics file: %s", output_csv)
        df = pd.read_csv(output_csv)
    else:
        all_rows = []
        for spec in [rpv_spec(root), dor_spec(root)]:
            drug_rows = drms[drms["drug"] == spec.structure.name.upper()]
            if drug_rows.empty:
                logging.warning("No DRM entries found for %s", spec.structure.name)
                continue
            chain_map = load_chain_subunits(spec.structure.cif_path)
            residue_maps = load_residue_mappings(spec.structure.cif_path)
            numbering = detect_numbering_scheme(spec.structure.cif_path, chain_map)
            all_rows.extend(
                run_mutations(
                    spec,
                    paths,
                    drug_rows,
                    chain_map,
                    residue_maps,
                    numbering,
                    replicates=args.replicates,
                    jitter_seed_base=args.seed,
                    jitter_angstrom=args.jitter_angstrom,
                    alchemy_config=alchemy_cfg,
                )
            )
        df = pd.DataFrame(all_rows)
    df.to_csv(output_csv, index=False)

    write_metrics_xlsx(df, paths.results / "metrics_summary.xlsx")
    plot_delta_metrics(df, paths)


if __name__ == "__main__":
    main()
