from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable
import re
import shlex


@dataclass(frozen=True)
class Paths:
    root: Path
    data: Path
    structures: Path
    ligands: Path
    generated: Path
    results: Path
    plots: Path


def project_paths(root: Path) -> Paths:
    data = root / "data"
    return Paths(
        root=root,
        data=data,
        structures=data / "structures",
        ligands=data / "ligands",
        generated=data / "generated",
        results=root / "results",
        plots=root / "results" / "plots",
    )


def ensure_dirs(paths: Iterable[Path]) -> None:
    for path in paths:
        path.mkdir(parents=True, exist_ok=True)


def sanitize_label(label: str) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "_", label).strip("_")


def parse_mutation_token(token: str) -> tuple[str, str, str]:
    match = re.match(r"^([A-Z])(\d+)([A-Z])$", token.strip().upper())
    if not match:
        raise ValueError(f"Unsupported mutation token: {token}")
    old_res, residue_id, new_res = match.groups()
    return old_res, residue_id, new_res


def parse_mutation_group(
    mutation: str, chains: str | list[str]
) -> list[tuple[str, str, str]]:
    tokens = mutation.split("+")
    if isinstance(chains, str):
        chain_list = [c.strip().upper() for c in chains.split("+") if c.strip()]
    else:
        chain_list = [str(c).strip().upper() for c in chains if str(c).strip()]

    if not chain_list:
        raise ValueError("No chains provided for mutation group.")
    if len(chain_list) == 1:
        chain_list = chain_list * len(tokens)
    if len(chain_list) != len(tokens):
        raise ValueError(
            f"Chain count ({len(chain_list)}) does not match mutation count ({len(tokens)})."
        )

    steps = []
    for token, chain_id in zip(tokens, chain_list):
        _, residue_id, new_res = parse_mutation_token(token)
        steps.append((chain_id, residue_id, new_res))
    return steps


def load_chain_subunits(cif_path: Path) -> dict[str, str]:
    lines = cif_path.read_text().splitlines()
    entity_names: dict[str, str] = {}
    chain_entities: dict[str, str] = {}

    idx = 0
    while idx < len(lines):
        line = lines[idx].strip()
        if line != "loop_":
            idx += 1
            continue

        idx += 1
        tags = []
        while idx < len(lines):
            tag_line = lines[idx].strip()
            if tag_line.startswith("_"):
                tags.append(tag_line)
                idx += 1
            else:
                break

        if not tags:
            continue

        data_tokens: list[str] = []
        while idx < len(lines):
            data_line = lines[idx]
            stripped = data_line.strip()
            if not stripped:
                idx += 1
                continue
            if stripped.startswith("#") or stripped == "loop_" or stripped.startswith("data_"):
                break
            if stripped.startswith("_"):
                break
            if stripped.startswith(";"):
                block_lines = [stripped[1:]]
                idx += 1
                while idx < len(lines):
                    block_line = lines[idx]
                    if block_line.startswith(";"):
                        block_lines.append(block_line[1:])
                        idx += 1
                        break
                    block_lines.append(block_line)
                    idx += 1
                data_tokens.append("\n".join(block_lines).strip())
                continue
            data_tokens.extend(shlex.split(stripped))
            idx += 1

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
