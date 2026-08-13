"""Canonical repository paths — the single source of truth for the on-disk layout.

Import these instead of hardcoding directory strings::

    from src.paths import MD_RUNS, FEP_PMX, FIGURES

Rationale: result/data paths are currently hardcoded across ~120 files, so any
move means editing all of them. Routing through this module makes the structure
explicit and turns a future move into a one-line change here.

Nothing in this module touches the filesystem on import; it only defines paths.
"""
from __future__ import annotations

from pathlib import Path

# src/paths.py  ->  parents[1] is the repository root.
REPO_ROOT = Path(__file__).resolve().parents[1]

# --- inputs (read-only provenance) -----------------------------------------
DATA = REPO_ROOT / "data"
STRUCTURES = DATA / "structures"          # experimental structures, e.g. 4NCG.cif
LIGANDS = DATA / "ligands"                # ligand definitions, e.g. dor.sdf
PREPARED = DATA / "prepared"              # prepared/parameterized inputs

# --- code & run provenance --------------------------------------------------
SRC = REPO_ROOT / "src"
SCRIPTS = REPO_ROOT / "scripts"
MANIFESTS = REPO_ROOT / "manifests"
LOGS = REPO_ROOT / "logs"

# --- results ----------------------------------------------------------------
RESULTS = REPO_ROOT / "results"
MD_RUNS = RESULTS / "md_runs"             # holo (DOR-bound) classical MD
APO_MD_RUNS = MD_RUNS / "apo"             # apo (ligand-free) classical MD
ANALYSIS = RESULTS / "analysis"
FEP_PMX = ANALYSIS / "fep_pmx"            # pmx non-equilibrium alchemical FEP

# --- human-facing outputs ---------------------------------------------------
FIGURES = REPO_ROOT / "figures"           # curated, manuscript-facing figures
MANUSCRIPT = REPO_ROOT / "manuscript"


def rel(path: Path | str) -> str:
    """Return `path` as a repo-relative POSIX string (for metadata/logging).

    Falls back to the absolute string if `path` is outside the repository.
    """
    p = Path(path).resolve()
    try:
        return p.relative_to(REPO_ROOT).as_posix()
    except ValueError:
        return p.as_posix()


def md_run_dir(genotype: str, replicate: int, *, apo: bool = False) -> Path:
    """Directory for one classical-MD run.

    Holo: results/md_runs/<genotype>/rep_NN
    Apo:  results/md_runs/apo/<genotype>/rep_NN
    """
    base = APO_MD_RUNS if apo else MD_RUNS
    return base / genotype / f"rep_{replicate:02d}"
