from __future__ import annotations

from pathlib import Path

from .cif_parser import iter_cif_loops


def load_residue_mappings(cif_path: Path) -> dict[str, dict[str, dict[str, str]]]:
    aa3 = {
        "ALA",
        "CYS",
        "ASP",
        "GLU",
        "PHE",
        "GLY",
        "HIS",
        "ILE",
        "LYS",
        "LEU",
        "MET",
        "ASN",
        "PRO",
        "GLN",
        "ARG",
        "SER",
        "THR",
        "VAL",
        "TRP",
        "TYR",
    }
    lines = cif_path.read_text().splitlines()
    maps: dict[str, dict[str, dict[str, str]]] = {}

    for tags, data_tokens in iter_cif_loops(lines):
        if not any(tag.startswith("_atom_site.") for tag in tags):
            continue

        def _idx(name: str) -> int | None:
            return tags.index(name) if name in tags else None

        auth_asym = _idx("_atom_site.auth_asym_id")
        auth_seq = _idx("_atom_site.auth_seq_id")
        auth_comp = _idx("_atom_site.auth_comp_id")
        label_asym = _idx("_atom_site.label_asym_id")
        label_seq = _idx("_atom_site.label_seq_id")
        label_comp = _idx("_atom_site.label_comp_id")

        ncols = len(tags)
        for row in range(0, len(data_tokens), ncols):
            values = data_tokens[row : row + ncols]
            if len(values) < ncols:
                break

            def _get(idx):
                return values[idx] if idx is not None else None

            auth_chain = _get(auth_asym)
            auth_seq_id = _get(auth_seq)
            auth_comp_id = _get(auth_comp)
            label_chain = _get(label_asym)
            label_seq_id = _get(label_seq)
            label_comp_id = _get(label_comp)

            if auth_chain and auth_seq_id and auth_comp_id:
                if auth_comp_id in aa3 and auth_seq_id not in {".", "?"}:
                    maps.setdefault(auth_chain, {}).setdefault("auth_map", {})[
                        auth_seq_id
                    ] = auth_comp_id
                if (
                    auth_comp_id in aa3
                    and auth_seq_id not in {".", "?"}
                    and label_seq_id
                    and label_seq_id not in {".", "?"}
                ):
                    maps.setdefault(auth_chain, {}).setdefault("auth_to_label", {})[
                        auth_seq_id
                    ] = label_seq_id

            if label_chain and label_seq_id and label_comp_id:
                if label_comp_id in aa3 and label_seq_id not in {".", "?"}:
                    maps.setdefault(label_chain, {}).setdefault("label_map", {})[
                        label_seq_id
                    ] = label_comp_id
                if (
                    label_comp_id in aa3
                    and label_seq_id not in {".", "?"}
                    and auth_seq_id
                    and auth_seq_id not in {".", "?"}
                ):
                    maps.setdefault(label_chain, {}).setdefault("label_to_auth", {})[
                        label_seq_id
                    ] = auth_seq_id

        break

    return maps
