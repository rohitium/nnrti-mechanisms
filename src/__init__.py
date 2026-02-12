"""NNRTI resistance analysis toolkit.

Canonical package layout:
- ``src.structure_prep``: mutation + preparation of WT/mutant systems
- ``src.md``: manifests, MD execution, OpenMM runtime, Sherlock helpers
- ``src.analysis``: susceptibility I/O, metrics, result collection, plotting
- ``src.utils``: shared filesystem/CIF/mutation utilities
"""

__all__ = [
    "analysis",
    "md",
    "structure_prep",
    "utils",
]
