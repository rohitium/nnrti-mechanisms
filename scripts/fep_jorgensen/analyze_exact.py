"""Analyze MCPRO mutation legs using the Jorgensen inhibitor-relative cycle."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path


def read_mutation_legs(path: Path) -> dict[str, tuple[float, float]]:
    """Read inhibitor,delta_g_mutation_kcal_mol,uncertainty_kcal_mol rows."""
    values: dict[str, tuple[float, float]] = {}
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        required = {
            "inhibitor", "delta_g_mutation_kcal_mol", "uncertainty_kcal_mol"
        }
        missing = required - set(reader.fieldnames or ())
        if missing:
            raise ValueError(f"Missing columns: {', '.join(sorted(missing))}")
        for row in reader:
            inhibitor = row["inhibitor"].strip().lower()
            if inhibitor in values:
                raise ValueError(f"Duplicate inhibitor: {inhibitor}")
            values[inhibitor] = (
                float(row["delta_g_mutation_kcal_mol"]),
                float(row["uncertainty_kcal_mol"]),
            )
    if "sustiva" not in values:
        raise ValueError("Exact Jorgensen normalization requires a sustiva leg")
    return values


def normalize_to_sustiva(
    legs: dict[str, tuple[float, float]], mutation: str
) -> list[dict[str, float | str]]:
    """Return ΔΔG_i = ΔG_mut(i) - ΔG_mut(Sustiva).

    Independent leg uncertainties are propagated in quadrature.  Sustiva is
    defined as exactly zero in the normalized table.
    """
    reference_dg, reference_sigma = legs["sustiva"]
    rows = []
    for inhibitor, (delta_g, sigma) in legs.items():
        is_reference = inhibitor == "sustiva"
        rows.append(
            {
                "mutation": mutation.upper(),
                "reference_inhibitor": "sustiva",
                "inhibitor": inhibitor,
                "delta_delta_g_kcal_mol": 0.0 if is_reference else delta_g - reference_dg,
                "uncertainty_kcal_mol": 0.0 if is_reference else math.hypot(sigma, reference_sigma),
            }
        )
    return rows


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Normalize MCPRO protein-mutation free energies to Sustiva"
    )
    parser.add_argument("input_csv", type=Path)
    parser.add_argument("--mutation", default="V106A")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    rows = normalize_to_sustiva(read_mutation_legs(args.input_csv), args.mutation)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with args.output.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
    print(json.dumps(rows, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
