"""Snakemake script: generate all publication plots."""

import logging
from pathlib import Path

import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

from src.analysis.plotting import (
    plot_all_metrics_vs_fold_reduction,
    plot_boundness_qc,
    plot_simulation_convergence,
)
from src.utils import project_paths

root = Path(".").resolve()
paths = project_paths(root)
paths.plots.mkdir(parents=True, exist_ok=True)

ddg_df = pd.read_csv(snakemake.input.ddg_full)  # noqa: F821
plot_all_metrics_vs_fold_reduction(ddg_df, paths)

boundness_df = pd.read_csv(snakemake.input.boundness)  # noqa: F821
plot_boundness_qc(boundness_df, paths)

rmsd_df = pd.read_csv(snakemake.input.rmsd)  # noqa: F821
com_df = pd.read_csv(snakemake.input.com)  # noqa: F821
plot_simulation_convergence(rmsd_df, com_df, paths)
