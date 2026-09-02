"""Molecular dynamics domain: manifest, worker, OpenMM helpers, and Sherlock tools."""

from .manifest import MDTask, get_task_by_id, load_manifest, save_manifest

__all__ = [
    "MDTask",
    "get_task_by_id",
    "load_manifest",
    "run_md_task",
    "save_manifest",
]


def __getattr__(name: str):
    if name == "run_md_task":
        from .worker import run_md_task

        return run_md_task
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
