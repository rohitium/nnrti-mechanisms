from __future__ import annotations

__all__ = [
    "MDTask",
    "collect_md_results",
    "compute_binding_ddg",
    "compute_correlations",
    "get_task_by_id",
    "load_manifest",
    "merge_with_structural_metrics",
    "run_result_collection",
    "save_manifest",
    "summarize_ddg_by_mutation",
]


def __getattr__(name: str):
    # Lazy imports keep optional dependencies (e.g. pandas) out of GPU worker startup.
    if name in {"MDTask", "get_task_by_id", "load_manifest", "save_manifest"}:
        from .manifest import MDTask, get_task_by_id, load_manifest, save_manifest

        return {
            "MDTask": MDTask,
            "get_task_by_id": get_task_by_id,
            "load_manifest": load_manifest,
            "save_manifest": save_manifest,
        }[name]

    if name in {
        "collect_md_results",
        "compute_binding_ddg",
        "compute_correlations",
        "merge_with_structural_metrics",
        "run_result_collection",
        "summarize_ddg_by_mutation",
    }:
        from .result_collector import (
            collect_md_results,
            compute_binding_ddg,
            compute_correlations,
            merge_with_structural_metrics,
            run_result_collection,
            summarize_ddg_by_mutation,
        )

        return {
            "collect_md_results": collect_md_results,
            "compute_binding_ddg": compute_binding_ddg,
            "compute_correlations": compute_correlations,
            "merge_with_structural_metrics": merge_with_structural_metrics,
            "run_result_collection": run_result_collection,
            "summarize_ddg_by_mutation": summarize_ddg_by_mutation,
        }[name]

    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
