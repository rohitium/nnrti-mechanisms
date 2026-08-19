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

# P1 single-residue legs for the ranking gate (Spearman vs experimental fold).
# Split by net charge: neutral legs run the standard pipeline now; the two
# charge-changing legs need the co-alchemical ion / double-box protocol (PLAN
# §6.2, not yet implemented) and are deferred.
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

# Co-alchemical ion ABANDONED 2026-08-05: decoupling one bulk Cl- to keep the box
# neutral dissipated ~20-26 kcal/mol (vs ~1-3 neutral), near-zero overlap, SEM ~1.4,
# BAR-Jarz disagreement ~3.7 — it does not converge (see OPERATIONS.md §7). We now
# run charge legs in a RAW non-neutral box (genion -neutral neutralises the A-state;
# the B-state carries net charge under PME's uniform background) and apply the
# Rocklin/Hunenberger analytical net-charge correction post-hoc (zero added
# perturbation). Set True only to reproduce the abandoned co-alchemical experiment.
USE_COALCHEMICAL_ION = False

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
# F227C (Phe→Cys) was TESTED at 500 ps (2026-08-08) then REVERTED to 100 ps
# (2026-08-10): the 500 ps run gave per-rep ΔG identical to 100 ps (holo
# −0.18/−1.72/−2.15 vs −0.1/−1.8/−2.2) — clean switch-length invariance. Its noise
# is NOT switch dissipation but across-rep endpoint scatter (fast local pocket
# repacking around the deleted ring; no slow basin seen in 100 ns MD once PBC
# artifacts are removed). The fix is more replicates (SEM∝1/√n), not longer
# switches, so F227C stays at the cheaper 100 ps default. The 3 existing 500 ps
# reps remain on disk and double as the switch-length-invariance check. Same
# reasoning applies to the other aromatic→Cys legs (V106I_to_V106I_F227C,
# wt_to_Y181C): keep them at 100 ps and add replicates.
# G190E was MOVED OUT of the long-switch set (2026-08-15). It was listed here
# speculatively (charge leg => assumed dissipative) back when charge legs still ran
# the co-alchemical ion. That ion is abandoned (OPERATIONS.md §7): the ~20-26 kcal
# dissipation came from annihilating a whole Cl-, not from the Gly->Glu mutation, and
# the raw-box + analytical-correction protocol has neutral-like dissipation. K103N —
# the other delta_q = -1 leg — ran the whole panel at 100 ps. G190E now matches it
# exactly, which also makes it a like-for-like member of the charge family for the
# SEM work. Switch length is in any case NOT the SEM lever (V106A 100 vs 500 ps:
# +1.69 +- 0.70 vs +1.76 +- 0.51, invariant); endpoint equilibration is.
_BASE_LONG_SWITCH_LEGS = {"wt_to_V106A", "wt_to_Y188L"}
_EXTRA_LONG_SWITCH_LEGS = {
    leg.strip() for leg in os.environ.get("NEQ_EXTRA_LONG_SWITCH_LEGS", "").split(",") if leg.strip()
}
LONG_SWITCH_LEGS = frozenset(_BASE_LONG_SWITCH_LEGS | _EXTRA_LONG_SWITCH_LEGS)
NEQ_SWITCH_PS_DEFAULT = _env_float("NEQ_SWITCH_PS_DEFAULT", 100.0)
LONG_SWITCH_PS = _env_float("NEQ_LONG_SWITCH_PS", 500.0)

# SLURM array bundling: snapshots executed sequentially per GPU job (see docs/pmx-neq-fep-plan.md §7).
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

