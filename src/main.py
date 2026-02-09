from __future__ import annotations

import argparse
import logging
from pathlib import Path

import pandas as pd

from .cluster import generate_slurm_script, run_result_collection
from .config import dor_4ncg_spec
from .dor_md_pipeline import prepare_local_openmm_only_for_cluster
from .openmm.md_protocol import MDProtocolConfig
from .plotting import (
    cleanup_legacy_plots,
    plot_all_metrics_vs_fold_reduction,
    plot_boundness_qc,
    plot_si_figure_s1_like,
    plot_simulation_convergence,
)
from .utils import ensure_dirs, project_paths


def main(argv: list[str] | None = None) -> None:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    logging.getLogger("MDAnalysis").setLevel(logging.WARNING)
    logging.getLogger("MDAnalysis.analysis").setLevel(logging.WARNING)

    parser = argparse.ArgumentParser(
        description="NNRTI explicit-MD workflow: Sherlock MD + local MM/GBSA analysis."
    )
    parser.add_argument("--prepare-local-openmm-only", action="store_true")
    parser.add_argument("--generate-slurm", action="store_true")
    parser.add_argument("--collect-results", action="store_true")

    parser.add_argument("--replicates", type=int, default=1)
    parser.add_argument("--jitter-angstrom", type=float, default=0.1)
    parser.add_argument("--seed", type=int, default=None)
    parser.add_argument("--mutation", type=str, default=None)

    parser.add_argument("--susceptibility-xlsx", type=Path, default=None)
    parser.add_argument("--prepared-dir", type=Path, default=None)
    parser.add_argument("--fep-manifest", type=Path, default=None)
    parser.add_argument("--manifest", type=Path, default=None)
    parser.add_argument("--slurm-script", type=Path, default=None)
    parser.add_argument("--fep-results-dir", type=Path, default=None)
    parser.add_argument("--results-dir", type=Path, default=None)

    parser.add_argument("--md-heating-ps", type=float, default=25.0)
    parser.add_argument("--md-production-ns", type=float, default=2.0)
    parser.add_argument("--trajectory-interval", type=int, default=2000)

    parser.add_argument("--slurm-partition", type=str, default="gpu")
    parser.add_argument("--slurm-time", type=str, default="6:00:00")
    parser.add_argument("--slurm-memory", type=str, default="16G")
    parser.add_argument("--conda-env", type=str, default=None)
    parser.add_argument("--use-openmm-module", action="store_true")

    parser.add_argument("--skip-structural-metrics", action="store_true")
    parser.add_argument("--metric-frame-stride", type=int, default=5)
    parser.add_argument("--metric-max-frames", type=int, default=200)
    parser.add_argument("--mmgbsa-snapshots", type=int, default=100)
    parser.add_argument("--mmgbsa-discard-fraction", type=float, default=0.25)

    args = parser.parse_args(argv)

    root = Path(__file__).resolve().parents[1]
    paths = project_paths(root)
    ensure_dirs([paths.generated, paths.results, paths.plots])

    susceptibility_xlsx = args.susceptibility_xlsx or (paths.data / "DRM-susceptibilities.csv.xlsx")
    prepared_dir = args.prepared_dir or (paths.data / "prepared" / "dor_4ncg")
    manifest = args.manifest or args.fep_manifest or (paths.results / "md_manifest.csv")
    slurm_script = args.slurm_script or (root / "scripts" / "sherlock" / "submit_all_tasks.sh")
    run_results_dir = args.results_dir or args.fep_results_dir or (paths.results / "md_runs")

    md_cfg = MDProtocolConfig(
        heating_ps=args.md_heating_ps,
        production_ns=args.md_production_ns,
        report_interval_steps=args.trajectory_interval,
    )

    if args.prepare_local_openmm_only:
        selected_mutations = {args.mutation} if args.mutation else None
        tasks = prepare_local_openmm_only_for_cluster(
            root=root,
            susceptibility_xlsx=susceptibility_xlsx,
            prepared_dir=prepared_dir,
            manifest_path=manifest,
            results_dir=run_results_dir,
            replicates=args.replicates,
            jitter_seed_base=args.seed,
            jitter_angstrom=args.jitter_angstrom,
            md_config=md_cfg,
            selected_mutations=selected_mutations,
        )
        n_structures = len({(t.mutation, t.replicate) for t in tasks})
        logging.info("Prepared %d structure/replicate setups and %d tasks", n_structures, len(tasks))
        logging.info("Manifest: %s", manifest)
        return

    if args.generate_slurm:
        if not manifest.exists():
            logging.error("Manifest not found: %s", manifest)
            return
        output = generate_slurm_script(
            manifest_path=manifest,
            output_script=slurm_script,
            partition=args.slurm_partition,
            time_limit=args.slurm_time,
            memory=args.slurm_memory,
            heating_ps=args.md_heating_ps,
            production_ns=args.md_production_ns,
            report_interval=args.trajectory_interval,
            conda_env=args.conda_env,
            use_openmm_module=args.use_openmm_module,
        )
        logging.info("Generated SLURM script: %s", output)
        return

    if args.collect_results:
        if not manifest.exists():
            logging.error("Manifest not found: %s", manifest)
            return

        cleanup_legacy_plots(paths)
        spec = dor_4ncg_spec(root)
        _, _, ddg_df = run_result_collection(
            manifest_path=manifest,
            fep_results_dir=run_results_dir,
            output_dir=paths.results,
            ligand_resname=spec.structure.ligand_resname,
            compute_structural=not args.skip_structural_metrics,
            metric_frame_stride=max(1, args.metric_frame_stride),
            metric_max_frames=max(1, args.metric_max_frames),
            mmgbsa_snapshots=max(5, args.mmgbsa_snapshots),
            mmgbsa_discard_fraction=max(0.0, min(0.9, args.mmgbsa_discard_fraction)),
        )
        logging.info("Wrote %s", paths.results / "ddg_summary.csv")
        logging.info("Wrote %s", paths.results / "correlation_analysis.csv")

        try:
            plot_all_metrics_vs_fold_reduction(ddg_df, paths)
            logging.info("Wrote %s", paths.plots / "all_metrics_vs_fold_reduction.png")
        except Exception as exc:
            logging.warning("Could not generate multi-metric plot: %s", exc)

        try:
            boundness_path = paths.results / "boundness_qc.csv"
            pos_summary_path = paths.results / "mutation_position_summary.csv"
            rmsd_path = paths.results / "rmsd_ca_profiles.csv"
            com_path = paths.results / "com_distance_profiles.csv"
            if boundness_path.exists():
                plot_boundness_qc(pd.read_csv(boundness_path), paths)
                logging.info("Wrote %s", paths.plots / "boundness_qc_min_distance.png")
            if pos_summary_path.exists():
                plot_si_figure_s1_like(pd.read_csv(pos_summary_path), paths)
                logging.info("Wrote %s", paths.plots / "fig_s1_like_mutation_landscape.png")
            rmsd_data = pd.read_csv(rmsd_path) if rmsd_path.exists() else pd.DataFrame()
            com_data = pd.read_csv(com_path) if com_path.exists() else pd.DataFrame()
            if not rmsd_data.empty or not com_data.empty:
                plot_simulation_convergence(rmsd_data, com_data, paths)
                logging.info("Wrote %s", paths.plots / "simulation_convergence.png")
        except Exception as exc:
            logging.warning("Could not generate selected plots: %s", exc)
        return

    parser.error("No action requested. Use one of: --prepare-local-openmm-only, --generate-slurm, --collect-results")


if __name__ == "__main__":
    main()
