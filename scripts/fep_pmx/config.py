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

# NEQ protocol (PLAN.md §4.3–4.4; pmx protein_mut tutorial)
NEQ_TEMPERATURE_K = 300.0
NEQ_DT_PS = 0.002
NEQ_EQUIL_NS = 5.0
NEQ_SWITCH_PS_DEFAULT = 100.0
# P0 panel: 12 systems × bundled switch tasks (36 at 100 snapshots; see docs/pmx-neq-fep-plan.md §7.2).
NEQ_SNAPSHOTS_DEFAULT = 100
NEQ_EQUIL_SNAPSHOT_START_PS = 100.0  # skip first 100 ps of equil when extracting
LONG_SWITCH_LEGS = frozenset({"wt_to_Y188L", "wt_to_G190E"})
LONG_SWITCH_PS = 500.0

# SLURM array bundling: snapshots executed sequentially per GPU job (see docs/pmx-neq-fep-plan.md §7).
SWITCH_SNAPSHOTS_PER_TASK_DEFAULT = 100
SWITCH_SNAPSHOTS_PER_TASK_LONG = 50  # 500 ps switches → ~15 h/task at 50 snaps


def switch_snapshots_per_task(leg_id: str) -> int:
    """Max NEQ switches run sequentially in one GPU array element."""
    if leg_id in LONG_SWITCH_LEGS:
        return SWITCH_SNAPSHOTS_PER_TASK_LONG
    return SWITCH_SNAPSHOTS_PER_TASK_DEFAULT


def switch_bundle_ranges(n_snapshots: int, leg_id: str) -> list[tuple[int, int]]:
    """Inclusive (start, end) snapshot index ranges for bundled switch tasks."""
    chunk = switch_snapshots_per_task(leg_id)
    ranges: list[tuple[int, int]] = []
    for start in range(0, n_snapshots, chunk):
        end = min(start + chunk, n_snapshots) - 1
        ranges.append((start, end))
    return ranges


def switch_ps_for_leg(leg_id: str) -> float:
    """Return NEQ switch length in ps for a leg."""
    if leg_id in LONG_SWITCH_LEGS:
        return LONG_SWITCH_PS
    return NEQ_SWITCH_PS_DEFAULT


def nsteps_for_time_ps(time_ps: float, *, dt_ps: float = NEQ_DT_PS) -> int:
    return max(1, int(round(time_ps / dt_ps)))


def delta_lambda_for_switch(switch_ps: float, *, dt_ps: float = NEQ_DT_PS) -> float:
    """Linear λ ramp 0→1 (or 1→0) over switch_ps."""
    return 1.0 / nsteps_for_time_ps(switch_ps, dt_ps=dt_ps)

