#!/usr/bin/env python3
"""Add body paragraphs after an existing paragraph, in its own formatting.

edit_docx_text.py rewrites text inside a paragraph that already exists. A
section that grows from one paragraph to several needs the other half: new
<w:p> elements, carrying the style of the paragraph they follow so the result
does not arrive in Word's default font next to the document's own.

The new paragraphs copy the <w:pPr> and the run properties of the anchor
paragraph's first run, and nothing else -- no fields, no bookmarks, no
comment anchors are duplicated.

Text is read from a UTF-8 file, one paragraph per blank-line-separated block,
so the source stays legible and reviewable next to the document.

Usage
-----
    python ops/docx/insert_paragraphs.py FILE.docx paragraphs.txt \
        --after "Three 100 ns apo replicates were run" [--out OUT.docx]
"""

from __future__ import annotations

import argparse
import copy
import re
import shutil
import zipfile
from pathlib import Path
import xml.etree.ElementTree as ET

W = "{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"
ET.register_namespace("w", W[1:-1])


def register_document_namespaces(xml: bytes) -> None:
    """Keep the document's own namespace prefixes across a round trip.

    ElementTree invents a prefix for every namespace it was not told about,
    so r:embed and wp:extent come back as ns3:embed and ns4:extent. Word reads
    the URI and does not care, but ops/docx/replace_figure.py matches on the
    prefix and stops finding the drawings in any document that has passed
    through here. Registering what the file declares keeps the edit local.
    """
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
    ap.add_argument("source", type=Path, help="UTF-8 text, blank line between paragraphs")
    ap.add_argument("--after", required=True,
                    help="leading text of the paragraph to insert after")
    ap.add_argument("--out", type=Path, default=None)
    args = ap.parse_args(argv)

    blocks = [b.strip() for b in args.source.read_text(encoding="utf8").split("\n\n")]
    blocks = [" ".join(b.split()) for b in blocks if b.strip()]
    if not blocks:
        raise SystemExit("no paragraphs in source")

    with zipfile.ZipFile(args.docx) as z:
        names = z.namelist()
        blobs = {n: z.read(n) for n in names}

    register_document_namespaces(blobs["word/document.xml"])
    root = ET.fromstring(blobs["word/document.xml"])
    body = root.find(f"{W}body")
    children = list(body)
    hits = [i for i, el in enumerate(children)
            if el.tag == W + "p" and _text(el).startswith(args.after)]
    if len(hits) != 1:
        raise SystemExit(f"REFUSED: {len(hits)} paragraphs start with "
                         f"{args.after[:60]!r}; it must match exactly one")
    at = hits[0]
    anchor = children[at]
    ppr = anchor.find(W + "pPr")
    first_run = anchor.find(W + "r")
    rpr = first_run.find(W + "rPr") if first_run is not None else None

    for offset, text in enumerate(blocks, start=1):
        par = ET.Element(W + "p")
        if ppr is not None:
            par.append(copy.deepcopy(ppr))
        run = ET.SubElement(par, W + "r")
        if rpr is not None:
            run.append(copy.deepcopy(rpr))
        t = ET.SubElement(run, W + "t")
        t.text = text
        t.set("{http://www.w3.org/XML/1998/namespace}space", "preserve")
        body.insert(at + offset, par)

    blobs["word/document.xml"] = ET.tostring(root, encoding="UTF-8",
                                             xml_declaration=True)
    out = args.out or args.docx
    if out == args.docx:
        shutil.copy2(args.docx, args.docx.with_suffix(args.docx.suffix + ".bak"))
    with zipfile.ZipFile(out, "w", zipfile.ZIP_DEFLATED) as z:
        for n in names:
            z.writestr(n, blobs[n])
    print(f"inserted {len(blocks)} paragraph(s) after {args.after[:50]!r} -> {out}")
    for b in blocks:
        print(f"   + {b[:70]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
