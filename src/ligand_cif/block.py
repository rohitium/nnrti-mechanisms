from __future__ import annotations

from pathlib import Path


def bond_order(value_order: str) -> int:
    mapping = {"SING": 1, "DOUB": 2, "TRIP": 3, "AROM": 4}
    return mapping.get(value_order.upper(), 1)


def normalize_element(element: str) -> str:
    element = element.strip()
    if not element:
        return element
    return element[0].upper() + element[1:].lower()


def load_block(cif_path: Path):
    import gemmi

    doc = gemmi.cif.read_file(str(cif_path))
    return doc.sole_block()
