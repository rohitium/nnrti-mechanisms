from __future__ import annotations

from pathlib import Path

from .cif_parser import iter_cif_loops


def load_chain_subunits(cif_path: Path) -> dict[str, str]:
    lines = cif_path.read_text().splitlines()
    entity_names: dict[str, str] = {}
    chain_entities: dict[str, str] = {}
    entity_types: dict[str, str] = {}

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
        elif "_entity_poly.entity_id" in tags and "_entity_poly.type" in tags:
            id_idx = tags.index("_entity_poly.entity_id")
            type_idx = tags.index("_entity_poly.type")
            ncols = len(tags)
            for row in range(0, len(data_tokens), ncols):
                values = data_tokens[row : row + ncols]
                if len(values) < ncols:
                    break
                entity_types[values[id_idx]] = values[type_idx]

    if not chain_entities:
        raise ValueError(f"Missing _struct_asym in {cif_path}")

    chain_map = {}
    for chain_id, entity_id in chain_entities.items():
        name = entity_names.get(str(entity_id), "")
        lower = name.lower()
        if "p66" in lower:
            chain_map[str(chain_id)] = "p66"
        elif "ribonuclease h" in lower:
            chain_map[str(chain_id)] = "p66"
        elif "p51" in lower or "p55" in lower:
            chain_map[str(chain_id)] = "p51"

    # Some structures (e.g. 4NCG) label p66 with a generic RT description.
    # If only one RT chain remains unmapped, infer it as p66.
    if any(v == "p51" for v in chain_map.values()):
        for chain_id, entity_id in chain_entities.items():
            if chain_id in chain_map:
                continue
            entity_type = entity_types.get(str(entity_id), "").lower()
            if "polypeptide" in entity_type:
                chain_map[str(chain_id)] = "p66"
                break
    return chain_map
