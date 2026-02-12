"""Analysis utilities: susceptibility IO, trajectory metrics, and plotting."""

from .metrics import ContactMetrics, EnsembleMetrics, compute_contacts, compute_ensemble_metrics, pocket_volume_proxy
from .plotting import (
    cleanup_legacy_plots,
    plot_all_metrics_vs_fold_reduction,
    plot_boundness_qc,
    plot_si_figure_s1_like,
    plot_simulation_convergence,
)
from .result_collector import (
    collect_ca_rmsd_profiles,
    collect_com_distance_profiles,
    collect_md_results,
    compute_binding_ddg,
    compute_boundness_qc,
    compute_correlations,
    compute_mmgbsa_metrics,
    compute_structural_metrics,
    merge_with_structural_metrics,
    run_result_collection,
)
from .susceptibility import load_dor_susceptibilities

__all__ = [
    "ContactMetrics",
    "EnsembleMetrics",
    "collect_ca_rmsd_profiles",
    "collect_com_distance_profiles",
    "collect_md_results",
    "cleanup_legacy_plots",
    "compute_binding_ddg",
    "compute_boundness_qc",
    "compute_contacts",
    "compute_correlations",
    "compute_ensemble_metrics",
    "compute_mmgbsa_metrics",
    "compute_structural_metrics",
    "load_dor_susceptibilities",
    "merge_with_structural_metrics",
    "plot_all_metrics_vs_fold_reduction",
    "plot_boundness_qc",
    "plot_si_figure_s1_like",
    "plot_simulation_convergence",
    "pocket_volume_proxy",
    "run_result_collection",
]
