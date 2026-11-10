#!/usr/bin/env python3
"""Remove whole paragraphs from a .docx, by the text they start with.

edit_docx_text.py rewrites text inside a paragraph and insert_paragraphs.py
adds new ones; this is the third case, a figure or a caption that has to
leave the document entirely.

Guard rails, since this removes content rather than changing it:

  * an anchor that does not match exactly one paragraph is refused, so a
    typo cannot delete the wrong caption or two of them;
  * ``--with-image`` also removes the paragraph immediately before the
    caption when that paragraph holds a drawing and no text, which is how a
    figure sits in this document -- but only then, so a caption preceded by
    body text leaves that text alone;
  * the media part itself is left in the package. Word ignores an unreferenced
    image, and deleting it would invalidate any other drawing that happens to
    point at the same part.

Usage
-----
    python ops/docx/delete_paragraphs.py FILE.docx \
        --starts-with "Figure 6:" --starts-with "Figure 7:" \
        --with-image [--out OUT.docx]
"""

from __future__ import annotations

import argparse
import re
import shutil
import zipfile
from pathlib import Path
import xml.etree.ElementTree as ET

W = "{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"
A = "{http://schemas.openxmlformats.org/drawingml/2006/main}"


def register_namespaces(xml: bytes) -> None:
    """Keep the document's own prefixes; see edit_docx_text.py for why."""
    # The whole file, not just the header: Word declares w14, wp14 and
    # friends on the drawing elements themselves, and a header-only scan
    # left those to be renamed ns5, ns6, ns8 on the way out.
    for prefix, uri in re.findall(rb'xmlns:([A-Za-z0-9]+)="([^"]+)"', xml):
        name = prefix.decode()
        # ElementTree reserves the ns<digits> form for the prefixes it
        # invents itself and raises on any attempt to register one. A
        # document that has already been through this tool carries them, so
        # skipping is what lets it be run twice.
        if re.fullmatch(r"ns\d+", name):
            continue
        ET.register_namespace(name, uri.decode())


def text_of(p) -> str:
    return "".join(t.text or "" for t in p.iter(W + "t"))


def has_drawing(p) -> bool:
    return p.find(f".//{A}blip") is not None


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("docx", type=Path)
    ap.add_argument("--starts-with", action="append", required=True,
                    help="leading text of a paragraph to delete; repeatable")
    ap.add_argument("--with-image", action="store_true",
                    help="also delete a preceding drawing-only paragraph")
    ap.add_argument("--out", type=Path)
    args = ap.parse_args(argv)

    out = args.out or args.docx
    if out != args.docx:
        shutil.copy(args.docx, out)

    with zipfile.ZipFile(args.docx) as z:
        parts = {n: z.read(n) for n in z.namelist()}
    register_namespaces(parts["word/document.xml"])
    root = ET.fromstring(parts["word/document.xml"])
    body = root.find(W + "body")
    kids = list(body)

    doomed = []
    for anchor in args.starts_with:
        hits = [i for i, p in enumerate(kids)
                if p.tag == W + "p" and text_of(p).strip().startswith(anchor)]
        if len(hits) != 1:
            raise SystemExit(f"refused: {anchor!r} matched {len(hits)} paragraphs")
        i = hits[0]
        doomed.append(i)
        if args.with_image and i > 0:
            prev = kids[i - 1]
            if prev.tag == W + "p" and has_drawing(prev) and not text_of(prev).strip():
                doomed.append(i - 1)
        print(f"deleting [{i}] {text_of(kids[i])[:70]!r}"
              + ("  + its image" if i - 1 in doomed else ""))

    for i in sorted(set(doomed), reverse=True):
        body.remove(kids[i])

    parts["word/document.xml"] = ET.tostring(root, xml_declaration=True,
                                             encoding="UTF-8")
    with zipfile.ZipFile(out, "w", zipfile.ZIP_DEFLATED) as z:
        for n, b in parts.items():
            z.writestr(n, b)
    print(f"removed {len(set(doomed))} paragraphs -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
