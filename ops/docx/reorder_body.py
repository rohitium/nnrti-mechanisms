#!/usr/bin/env python3
"""Reorder the top-level blocks of a .docx body to a given sequence.

Restructuring the Supporting Information means moving sections and figures
without rebuilding the file, because rebuilding would lose the embedded
images. This splits the body into its top-level blocks -- paragraphs and
tables, in document order -- lets a plan reorder them by matching text, and
reassembles. Every block is preserved: a block the plan does not mention
keeps its position relative to the block it follows, and the count is checked
before anything is written.

The plan is a text file, one entry per line:
    KEEP <text>     the block whose text contains <text>, and any blocks
                    attached to it (a table immediately following a paragraph)
    DROP <text>     remove that block
    --              blank separator, ignored

Usage
-----
    python ops/docx/reorder_body.py FILE.docx --plan plan.txt [--dry-run]
"""

from __future__ import annotations

import argparse
import html
import os
import re
import sys
import zipfile
from pathlib import Path

DOC = "word/document.xml"
BLOCK_RE = re.compile(r"<w:tbl>.*?</w:tbl>|<w:p\b[^>]*/>|<w:p\b[^>]*>.*?</w:p>", re.S)
TEXT_RE = re.compile(r"<w:t(?: [^>]*)?>(.*?)</w:t>", re.S)


def block_text(b: str) -> str:
    b = re.sub(r"<w:instrText[^>]*>.*?</w:instrText>", "", b, flags=re.S)
    return re.sub(r"\s+", " ", html.unescape("".join(TEXT_RE.findall(b)))).strip()


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("docx", type=Path)
    ap.add_argument("--plan", type=Path, required=True)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args(argv)

    zin = zipfile.ZipFile(args.docx)
    items = {n: zin.read(n) for n in zin.namelist()}
    infos = zin.infolist()
    zin.close()
    doc = items[DOC].decode("utf-8")

    m = re.search(r"<w:body>(.*)</w:body>", doc, re.S)
    if not m:
        raise SystemExit("no <w:body> found")
    body = m.group(1)
    # trailing section properties must stay last
    sect = ""
    s = re.search(r"(<w:sectPr\b.*?</w:sectPr>)\s*$", body, re.S)
    if s:
        sect = s.group(1)
        body = body[:s.start()]

    blocks = BLOCK_RE.findall(body)
    if "".join(blocks) != body.strip():
        leftover = len(body.strip()) - len("".join(blocks))
        raise SystemExit(f"block split is lossy ({leftover} chars unaccounted); "
                         "refusing to reorder")
    print(f"{len(blocks)} top-level blocks")

    plan = [l.strip() for l in args.plan.read_text().splitlines()
            if l.strip() and not l.strip().startswith("--")]

    used, ordered, dropped = set(), [], []
    for line in plan:
        op, _, needle = line.partition(" ")
        needle = needle.strip()
        hits = [i for i, b in enumerate(blocks)
                if i not in used and needle and needle in block_text(b)]
        if not hits:
            raise SystemExit(f"no unused block matches {needle!r}")
        i = hits[0]
        used.add(i)
        group = [blocks[i]]
        # Blocks that follow with no text of their own belong to this one: a
        # table, the paragraph holding a figure's image, and the empty
        # paragraphs that space a figure out. Without this an image-only
        # paragraph is left behind and the figure is orphaned.
        j = i + 1
        while j < len(blocks) and j not in used and (
                blocks[j].startswith("<w:tbl>") or not block_text(blocks[j])):
            used.add(j); group.append(blocks[j]); j += 1
        if op == "DROP":
            dropped.extend(group)
        else:
            ordered.extend(group)

    tail = [b for i, b in enumerate(blocks) if i not in used]
    new_blocks = ordered + tail
    print(f"  placed {len(ordered)}, dropped {len(dropped)}, "
          f"kept in place {len(tail)}  -> {len(new_blocks)}")
    if len(new_blocks) + len(dropped) != len(blocks):
        raise SystemExit("block accounting failed; refusing to write")

    if args.dry_run:
        for b in new_blocks[:40]:
            t = block_text(b)
            print(f"   {'[TBL] ' if b.startswith('<w:tbl>') else ''}{t[:88]}")
        return 0

    new_body = "".join(new_blocks) + sect
    doc = doc[:m.start(1)] + new_body + doc[m.end(1):]
    items[DOC] = doc.encode("utf-8")
    tmp = args.docx.with_suffix(".tmp")
    with zipfile.ZipFile(tmp, "w", zipfile.ZIP_DEFLATED) as z:
        for info in infos:
            z.writestr(info, items[info.filename])
    os.replace(tmp, args.docx)
    print(f"wrote {args.docx}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
