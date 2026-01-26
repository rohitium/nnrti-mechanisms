from __future__ import annotations

from pathlib import Path

from .cif_parser import iter_cif_loops


def load_chain_subunits(cif_path: Path) -> dict[str, str]:
    lines = cif_path.read_text().splitlines()
    entity_names: dict[str, str] = {}
    chain_entities: dict[str, str] = {}

    for tags, data_tokens in iter_cif_loops(lines):
        if "_entity_name_com.entity_id" in tags and "_entity_name_com.name" in tags:
            id_idx = tags.index("_entity_name_com.entity_id")
            name_idx = tags.index("_entity_name_com.name")
            ncols = len(tags)
            for row in range(0, len(data_tokens), ncols):
                values = data_tokens[row : row + ncols]
                if len(values) < ncols:
                    break
                entity_names[values[id_idx]] = values[name_idx]
        elif "_entity.id" in tags and "_entity.pdbx_description" in tags:
            id_idx = tags.index("_entity.id")
            name_idx = tags.index("_entity.pdbx_description")
            ncols = len(tags)
            for row in range(0, len(data_tokens), ncols):
                values = data_tokens[row : row + ncols]
                if len(values) < ncols:
                    break
                entity_names[values[id_idx]] = values[name_idx]
        elif "_struct_asym.id" in tags and "_struct_asym.entity_id" in tags:
            id_idx = tags.index("_struct_asym.id")
            entity_idx = tags.index("_struct_asym.entity_id")
            ncols = len(tags)
            for row in range(0, len(data_tokens), ncols):
                values = data_tokens[row : row + ncols]
                if len(values) < ncols:
                    break
                chain_entities[values[id_idx]] = values[entity_idx]

    if not chain_entities:
        raise ValueError(f"Missing _struct_asym in {cif_path}")

    chain_map = {}
    for chain_id, entity_id in chain_entities.items():
        name = entity_names.get(str(entity_id), "")
        lower = name.lower()
        if "p66" in lower:
            chain_map[str(chain_id)] = "p66"
        elif "p51" in lower or "p55" in lower:
            chain_map[str(chain_id)] = "p51"
    return chain_map
