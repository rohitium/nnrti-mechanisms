from __future__ import annotations

import argparse
import logging
from pathlib import Path

import pandas as pd

from .config import dor_spec, rpv_spec
from .drm_io import load_drms
from .metrics_io import write_metrics_xlsx
from .mutation.runner import run_mutations
from .numbering import detect_numbering_scheme
from .plotting import plot_delta_metrics
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
    args = parser.parse_args(argv)

    root = Path(__file__).resolve().parents[1]
    paths = project_paths(root)
    ensure_dirs([paths.generated, paths.results, paths.plots])

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
                )
            )
        df = pd.DataFrame(all_rows)
    df.to_csv(output_csv, index=False)

    write_metrics_xlsx(df, paths.results / "metrics_summary.xlsx")
    plot_delta_metrics(df, paths)


if __name__ == "__main__":
    main()
