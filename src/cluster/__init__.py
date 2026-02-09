from __future__ import annotations

__all__ = [
    "MDTask",
    "FEPTask",
    "collect_fep_results",
    "compute_binding_ddg",
    "compute_correlations",
    "generate_slurm_script",
    "get_task_by_id",
    "get_task_count",
    "load_manifest",
    "merge_with_structural_metrics",
    "run_result_collection",
    "save_manifest",
    "summarize_ddg_by_mutation",
]


def __getattr__(name: str):
    # Lazy imports keep optional dependencies (e.g. pandas) out of GPU worker startup.
    if name in {"MDTask", "FEPTask", "get_task_by_id", "load_manifest", "save_manifest"}:
        from .manifest import MDTask, FEPTask, get_task_by_id, load_manifest, save_manifest

        return {
            "MDTask": MDTask,
            "FEPTask": FEPTask,
            "get_task_by_id": get_task_by_id,
            "load_manifest": load_manifest,
            "save_manifest": save_manifest,
        }[name]

    if name in {"generate_slurm_script", "get_task_count"}:
        from .slurm_generator import generate_slurm_script, get_task_count

        return {
            "generate_slurm_script": generate_slurm_script,
            "get_task_count": get_task_count,
        }[name]

    if name in {
        "collect_fep_results",
        "compute_binding_ddg",
        "compute_correlations",
        "merge_with_structural_metrics",
        "run_result_collection",
        "summarize_ddg_by_mutation",
    }:
        from .result_collector import (
            collect_fep_results,
            compute_binding_ddg,
            compute_correlations,
            merge_with_structural_metrics,
            run_result_collection,
            summarize_ddg_by_mutation,
        )

        return {
            "collect_fep_results": collect_fep_results,
            "compute_binding_ddg": compute_binding_ddg,
            "compute_correlations": compute_correlations,
            "merge_with_structural_metrics": merge_with_structural_metrics,
            "run_result_collection": run_result_collection,
            "summarize_ddg_by_mutation": summarize_ddg_by_mutation,
        }[name]

    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
