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
# C-rescale barostat warmup before Parrinello-Rahman production. Starting P-R
# directly from a minimized structure with generated velocities can blow up a
# large solvated box; a short C-rescale phase relaxes the box first.
NEQ_WARMUP_PS = 500.0
# Endpoint equilibration (ns), applied to BOTH endpoints. Env-overridable so a
# longer-equilibration sensitivity test needs no code edit — the P0 pilot showed
# the reverse (λ=1 / mutant) work distributions are the noisy half, so more
# equilibration is the second lever after switch length (docs/pmx-neq-fep-plan.md
# §3.4). Example: NEQ_EQUIL_NS=10 REPLICATES=3 FORCE=1 bash prepare_p0_neq.sh
NEQ_EQUIL_NS = _env_float("NEQ_EQUIL_NS", 5.0)
# Snapshots per endpoint per replicate. Increasing this tightens the statistical
# error of ΔG but does NOT improve Crooks overlap (same distributions) — use
# switch length / equilibration for overlap, snapshots for error bars.
NEQ_SNAPSHOTS_DEFAULT = 100
NEQ_EQUIL_SNAPSHOT_START_PS = 100.0  # skip first 100 ps of equil when extracting

# Legs whose driven switches run LONG_SWITCH_PS instead of NEQ_SWITCH_PS_DEFAULT.
# V106A was added after the P0 pilot: its 100 ps switches gave marginal Crooks
# overlap (large reverse-work dissipation, ~1.5–6 kcal fwd/rev gap); 500 ps
# reduces dissipation and improves overlap. Extend the set for further tests via
# NEQ_EXTRA_LONG_SWITCH_LEGS (comma-separated leg ids) without editing code.
_BASE_LONG_SWITCH_LEGS = {"wt_to_V106A", "wt_to_Y188L", "wt_to_G190E"}
_EXTRA_LONG_SWITCH_LEGS = {
    leg.strip() for leg in os.environ.get("NEQ_EXTRA_LONG_SWITCH_LEGS", "").split(",") if leg.strip()
}
LONG_SWITCH_LEGS = frozenset(_BASE_LONG_SWITCH_LEGS | _EXTRA_LONG_SWITCH_LEGS)
NEQ_SWITCH_PS_DEFAULT = _env_float("NEQ_SWITCH_PS_DEFAULT", 100.0)
LONG_SWITCH_PS = _env_float("NEQ_LONG_SWITCH_PS", 500.0)

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

