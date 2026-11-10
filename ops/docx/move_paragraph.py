#!/usr/bin/env python3
"""Move one paragraph to sit immediately before another.

Inserting new figures into a section pushes the figure that was already there
down past them, so Figure 5 ended up printed after Figures 6 and 7. Rebuilding
the body with reorder_body.py needs a plan naming every block; this moves a
single block and leaves everything else where it is.

The paragraph is moved whole, so a figure whose drawing and caption share one
<w:p> -- which is how the manuscript lays them out -- travels with its art.

Usage
-----
    python ops/docx/move_paragraph.py FILE.docx \
        --move "Figure 5: Mutation-induced" --before "Figure 6: RT-specific"
"""
from __future__ import annotations

import argparse
import re
import shutil
import zipfile
from pathlib import Path
import xml.etree.ElementTree as ET

W = "{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"
ET.register_namespace("w", W[1:-1])


def register_document_namespaces(xml: bytes) -> None:
    for prefix, uri in re.findall(
            r'xmlns:([A-Za-z0-9_.-]+)="([^"]+)"', xml.decode("utf8", "replace")):
        try:
            ET.register_namespace(prefix, uri)
        except ValueError:
            pass


def _text(par) -> str:
    return "".join(t.text or "" for t in par.iter(W + "t"))


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("docx", type=Path)
    ap.add_argument("--move", required=True, help="leading text of the paragraph to move")
    ap.add_argument("--before", required=True, help="leading text of its new successor")
    ap.add_argument("--out", type=Path, default=None)
    args = ap.parse_args(argv)

    with zipfile.ZipFile(args.docx) as z:
        names = z.namelist()
        blobs = {n: z.read(n) for n in names}
    register_document_namespaces(blobs["word/document.xml"])
    root = ET.fromstring(blobs["word/document.xml"])
    body = root.find(f"{W}body")
    kids = list(body)

    def find(prefix):
        hits = [i for i, el in enumerate(kids)
                if el.tag == W + "p" and _text(el).strip().startswith(prefix)]
        if len(hits) != 1:
            raise SystemExit(f"REFUSED: {len(hits)} paragraphs start with {prefix!r}")
        return hits[0]

    src, dst = find(args.move), find(args.before)
    el = kids[src]
    body.remove(el)
    # Recompute after the removal: indices above src have shifted down by one.
    kids2 = list(body)
    dst2 = [i for i, e in enumerate(kids2)
            if e.tag == W + "p" and _text(e).strip().startswith(args.before)][0]
    body.insert(dst2, el)

    blobs["word/document.xml"] = ET.tostring(root, encoding="UTF-8",
                                             xml_declaration=True)
    out = args.out or args.docx
    if out == args.docx:
        shutil.copy2(args.docx, args.docx.with_suffix(args.docx.suffix + ".bak"))
    with zipfile.ZipFile(out, "w", zipfile.ZIP_DEFLATED) as z:
        for n in names:
            z.writestr(n, blobs[n])
    print(f"moved {args.move!r} to before {args.before!r} -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
