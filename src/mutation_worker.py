"""Backward-compatible import for multiprocessing pickling."""

from .mutation.worker import mutation_worker

__all__ = ["mutation_worker"]
