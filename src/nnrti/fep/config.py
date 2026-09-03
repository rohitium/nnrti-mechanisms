"""Shared configuration for pmx + GROMACS NEQ FEP."""

from __future__ import annotations

import os
from pathlib import Path


def _env_float(name: str, default: float) -> float:
    raw = os.environ.get(name)
    if raw is None:
        return default
    try:
        return float(raw)
    except ValueError:
        return default


def _env_int(name: str, default: int) -> int:
    raw = os.environ.get(name)
    if raw is None:
        return default
    try:
        value = int(raw)
    except ValueError:
        return default
    return value if value > 0 else default


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

# Single-residue legs, split by whether the substitution changes net charge.
P1_NEUTRAL_LEGS = (
    "wt_to_F227C", "wt_to_G190A", "wt_to_V106I",
    "wt_to_V106M", "wt_to_Y181C", "wt_to_Y318F",
)
P1_CHARGE_LEGS = ("wt_to_K103N", "wt_to_G190E")

# Charge-changing legs: protein net-charge change (A-state -> B-state). Only the
# WT->charged single legs change net charge — the compound second legs (M230L etc.)
# are neutral.
CHARGE_LEG_DELTA_Q = {
    "wt_to_K103N": -1,   # Lys+ -> Asn0
    "wt_to_G190E": -1,   # Gly0 -> Glu-
}

# Charge-changing legs run in a non-neutral box: genion neutralises the A state,
# the B state carries net charge under PME's uniform background, and the
# Rocklin/Hunenberger analytical net-charge correction is applied afterwards.

# Match OpenMM MD prep (md_protocol.py / manuscript)
SOLVENT_PADDING_NM = 1.0
IONIC_STRENGTH_M = 0.15
BOX_TYPE = "dodecahedron"
WATER_MODEL = "tip3p"

# NEQ protocol (pmx protein_mut tutorial)
NEQ_TEMPERATURE_K = 300.0
NEQ_DT_PS = 0.002
# C-rescale barostat warmup before Parrinello-Rahman production. Starting P-R
# directly from a minimized structure with generated velocities can blow up a
# large solvated box; a short C-rescale phase relaxes the box first.
NEQ_WARMUP_PS = 500.0
# Endpoint equilibration (ns), applied to BOTH endpoints. Env-overridable so a
# Set NEQ_EQUIL_NS to lengthen end-state equilibration without editing code.
NEQ_EQUIL_NS = _env_float("NEQ_EQUIL_NS", 5.0)
# Snapshots per endpoint per replicate. Increasing this tightens the statistical
# error of ΔG but does NOT improve Crooks overlap (same distributions) — use
# switch length / equilibration for overlap, snapshots for error bars.
NEQ_SNAPSHOTS_DEFAULT = 100
NEQ_EQUIL_SNAPSHOT_START_PS = 100.0  # skip first 100 ps of equil when extracting

# Legs whose driven switches run LONG_SWITCH_PS instead of NEQ_SWITCH_PS_DEFAULT:
# those whose forward and reverse work distributions are widely separated at the
# default length. Extend the set via NEQ_EXTRA_LONG_SWITCH_LEGS (comma-separated
# leg ids) without editing code.
NEQ_SWITCH_PS_DEFAULT = _env_float("NEQ_SWITCH_PS_DEFAULT", 100.0)
LONG_SWITCH_PS = _env_float("NEQ_LONG_SWITCH_PS", 500.0)

# SLURM array bundling: snapshots executed sequentially per GPU job.
#
# Switches are embarrassingly parallel -- each is an independent run from its own
# snapshot -- so this constant, not GPU availability, sets the wall clock. A leg is
# 1200 switches; at 50/task that is 24 elements and ~15 h each. Lowering it trades
# more array elements (against the MaxSubmitPU=100 cap) for a shorter critical path.
#
# Both are env-overridable so a tail-end top-up can be re-bundled small and finish in
# hours instead of a day: on 2026-08-17 G190E's last 96 switches sat in 2 elements at
# ~21 h; re-bundled at 6/task they spread over 16 GPUs (~3 h). Re-bundling is SAFE --
# run_neq_task skips per-switch on an existing dgdl.xvg, not per-task, so completed
# work is never redone. Rebuild the manifest with prepare_neq --force afterwards.
SWITCH_SNAPSHOTS_PER_TASK_DEFAULT = _env_int("NEQ_SNAPSHOTS_PER_TASK", 100)
SWITCH_SNAPSHOTS_PER_TASK_LONG = _env_int("NEQ_SNAPSHOTS_PER_TASK_LONG", 50)


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

