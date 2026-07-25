"""Shared configuration for pmx + GROMACS NEQ FEP."""

from __future__ import annotations

from pathlib import Path

FEP_PMX_ROOT = Path("results/analysis/fep_pmx")
LEG_INPUTS = FEP_PMX_ROOT / "inputs"
DOR_ITP_DIR = LEG_INPUTS / "dor_openff"

PMX_FORCE_FIELD = "amber14sbmut"
# pdb2gmx -ff must match the *.ff directory name under GMXLIB/mutff (not amber14SB label).
PMX_FORCE_FIELD_LABEL = "amber14sbmut"
OPENFF_FORCE_FIELD = "openff-2.0.0.offxml"
LIGAND_SDF = Path("data/ligands/dor.sdf")
LIGAND_RESNAME = "2KW"

P0_LEGS = ("wt_to_V106A", "wt_to_Y188L")

# Match OpenMM MD prep (md_protocol.py / manuscript)
SOLVENT_PADDING_NM = 1.0
IONIC_STRENGTH_M = 0.15
BOX_TYPE = "dodecahedron"
WATER_MODEL = "tip3p"
