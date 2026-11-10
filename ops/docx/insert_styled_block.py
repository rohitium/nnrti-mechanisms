#!/usr/bin/env python3
"""Insert several paragraphs after an anchor, each in a style of its own.

insert_paragraphs.py gives every new paragraph the anchor's formatting, which
is right when a section grows by a paragraph and wrong when the new material
is a section: a heading, body text and a figure caption are three different
styles, and a heading set in caption type looks like a mistake.

Here each block names the paragraph it should be formatted like, by that
paragraph's leading text, so the new section reuses styles the document
already defines instead of inventing any.

The spec is JSON: [{"like": "9. Conformation", "text": "10. ..."}, ...],
applied in order after the anchor.

Usage
-----
    python ops/docx/insert_styled_block.py FILE.docx spec.json \
        --after "Supplementary Figure 6." [--out OUT.docx]
"""

from __future__ import annotations

import argparse
import copy
import json
import re
import shutil
import zipfile
from pathlib import Path
import xml.etree.ElementTree as ET

W = "{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"


def register_namespaces(xml: bytes) -> None:
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


def find_one(kids, anchor):
    hits = [p for p in kids
            if p.tag == W + "p" and text_of(p).strip().startswith(anchor)]
    if len(hits) != 1:
        raise SystemExit(f"refused: {anchor!r} matched {len(hits)} paragraphs")
    return hits[0]


def build(model, text):
    """A new paragraph carrying the model's pPr and first run's rPr."""
    p = ET.Element(W + "p")
    ppr = model.find(W + "pPr")
    if ppr is not None:
        p.append(copy.deepcopy(ppr))
    r = ET.SubElement(p, W + "r")
    model_r = model.find(W + "r")
    if model_r is not None:
        rpr = model_r.find(W + "rPr")
        if rpr is not None:
            r.append(copy.deepcopy(rpr))
    t = ET.SubElement(r, W + "t")
    t.set("{http://www.w3.org/XML/1998/namespace}space", "preserve")
    t.text = text
    return p


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("docx", type=Path)
    ap.add_argument("spec", type=Path)
    ap.add_argument("--after", required=True)
    ap.add_argument("--out", type=Path)
    args = ap.parse_args(argv)

    out = args.out or args.docx
    with zipfile.ZipFile(args.docx) as z:
        parts = {n: z.read(n) for n in z.namelist()}
    register_namespaces(parts["word/document.xml"])
    root = ET.fromstring(parts["word/document.xml"])
    body = root.find(W + "body")
    kids = list(body)

    anchor = find_one(kids, args.after)
    at = list(body).index(anchor) + 1
    for block in json.loads(args.spec.read_text()):
        model = find_one(kids, block["like"])
        body.insert(at, build(model, block["text"]))
        print(f"+ [{at}] like {block['like'][:28]!r}: {block['text'][:60]}")
        at += 1

    parts["word/document.xml"] = ET.tostring(root, xml_declaration=True,
                                             encoding="UTF-8")
    with zipfile.ZipFile(out, "w", zipfile.ZIP_DEFLATED) as z:
        for n, b in parts.items():
            z.writestr(n, b)
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
