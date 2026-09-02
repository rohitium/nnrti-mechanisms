from __future__ import annotations

import re
import zlib


def deterministic_seed(jitter_seed_base: int | None, safe_label: str, replicate: int) -> int | None:
    """Derive a per-replicate seed from a base seed, mutation label, and replicate number."""
    if jitter_seed_base is None:
        return None
    token = f"{safe_label}:{replicate}".encode()
    return int(jitter_seed_base + (zlib.crc32(token) % 100000))


def sanitize_label(label: str) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "_", label).strip("_")


def parse_mutation_token(token: str) -> tuple[str, str, str]:
    match = re.match(r"^([A-Z])(\d+)([A-Z])$", token.strip().upper())
    if not match:
        raise ValueError(f"Unsupported mutation token: {token}")
    old_res, residue_id, new_res = match.groups()
    return old_res, residue_id, new_res


def one_to_three(res_name: str) -> str:
    res_name = res_name.strip().upper()
    if len(res_name) == 3:
        return res_name
    one_to_three_map = {
        "A": "ALA",
        "C": "CYS",
        "D": "ASP",
        "E": "GLU",
        "F": "PHE",
        "G": "GLY",
        "H": "HIS",
        "I": "ILE",
        "K": "LYS",
        "L": "LEU",
        "M": "MET",
        "N": "ASN",
        "P": "PRO",
        "Q": "GLN",
        "R": "ARG",
        "S": "SER",
        "T": "THR",
        "V": "VAL",
        "W": "TRP",
        "Y": "TYR",
    }
    if res_name not in one_to_three_map:
        raise ValueError(f"Unknown residue code: {res_name}")
    return one_to_three_map[res_name]
