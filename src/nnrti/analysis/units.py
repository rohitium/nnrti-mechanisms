"""Energy-unit handling for the analysis pipeline.

Manuscript-facing energies are reported in **kcal/mol**, matching the pmx
non-equilibrium FEP outputs in ``scripts/fep_pmx`` (see ``analyze_neq.py``,
which uses the same ``4.184`` factor and an ``OUTPUT_UNITS`` marker).

OpenMM reports in kJ/mol, so ``src/nnrti/md/openmm/mmgbsa.py`` emits kJ/mol. The
conversion happens once, at the canonical-table rebuild boundary, and the
result is stamped into an ``energy_units`` column so the conversion is
idempotent and every consumer can tell what it is reading.
"""

from __future__ import annotations

import pandas as pd

KJ_PER_KCAL = 4.184
KJ_TO_KCAL = 1.0 / KJ_PER_KCAL

ENERGY_UNITS_COLUMN = "energy_units"
KCAL_UNITS = "kcal/mol"
KJ_UNITS = "kJ/mol"

#: MM/GBSA replicate-level energy columns, including their std/sem partners.
MMGBSA_ENERGY_COLUMNS: tuple[str, ...] = tuple(
    f"{base}{suffix}"
    for base in (
        "binding_dg",
        "binding_dg_vdw",
        "binding_dg_electrostatic",
        "binding_dg_gb",
        "binding_dg_sa",
    )
    for suffix in ("", "_std", "_sem")
)


def frame_energy_units(df: pd.DataFrame) -> str:
    """Units currently carried by ``df``; assumes kJ/mol when unstamped."""
    if ENERGY_UNITS_COLUMN not in df.columns:
        return KJ_UNITS
    values = df[ENERGY_UNITS_COLUMN].dropna().unique()
    if len(values) == 0:
        return KJ_UNITS
    if len(values) > 1:
        raise ValueError(f"Mixed {ENERGY_UNITS_COLUMN} values in frame: {sorted(values)}")
    return str(values[0])


def convert_energy_columns(df: pd.DataFrame, target_units: str, columns=None) -> pd.DataFrame:
    """Convert energy columns to ``target_units``, stamping the units column.

    Idempotent: a frame already stamped with ``target_units`` is returned with
    its values untouched, so re-running a rebuild cannot double-convert.
    """
    if target_units not in {KCAL_UNITS, KJ_UNITS}:
        raise ValueError(f"target_units must be {KCAL_UNITS!r} or {KJ_UNITS!r}, got {target_units!r}")

    out = df.copy()
    current = frame_energy_units(out)
    if current == target_units:
        out[ENERGY_UNITS_COLUMN] = target_units
        return out

    factor = KJ_TO_KCAL if target_units == KCAL_UNITS else KJ_PER_KCAL
    if columns is None:
        columns = MMGBSA_ENERGY_COLUMNS
    for column in columns:
        if column in out.columns:
            out[column] = pd.to_numeric(out[column], errors="coerce") * factor
    out[ENERGY_UNITS_COLUMN] = target_units
    return out


def read_energy_table(path, target_units: str = KCAL_UNITS, columns=None) -> pd.DataFrame:
    """Read an energy CSV and return it in ``target_units``.

    Frames written before the units stamp existed are assumed to be kJ/mol
    (the OpenMM-native output of ``src/nnrti/md/openmm/mmgbsa.py``) and are converted;
    frames already stamped with ``target_units`` pass through untouched.
    """
    return convert_energy_columns(pd.read_csv(path), target_units, columns=columns)
