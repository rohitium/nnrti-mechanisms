from __future__ import annotations

from typing import Iterator
import shlex


def iter_cif_loops(lines: list[str]) -> Iterator[tuple[list[str], list[str]]]:
    idx = 0
    while idx < len(lines):
        line = lines[idx].strip()
        if line != "loop_":
            idx += 1
            continue
        idx += 1
        tags: list[str] = []
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
            if stripped.startswith("#") or stripped == "loop_" or stripped.startswith(
                "data_"
            ):
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

        yield tags, data_tokens
