"""NNRTI resistance analysis toolkit.

Canonical package layout:
- ``nnrti.structure_prep``: mutation + preparation of WT/mutant systems
- ``nnrti.md``: manifests, MD execution, OpenMM runtime, Sherlock helpers
- ``nnrti.analysis``: susceptibility I/O, metrics, result collection, plotting
- ``nnrti.utils``: shared filesystem/CIF/mutation utilities
"""

__all__ = [
    "analysis",
    "md",
    "structure_prep",
    "utils",
]
