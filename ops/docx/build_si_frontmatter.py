#!/usr/bin/env python3
"""Give the Supporting Information a title page, a Contents page and headings.

JCIM Supporting Information opens with a header page carrying the manuscript
title and authors, then a Contents listing of the sections, then the numbered
sections themselves. The sections here were single paragraphs whose first
sentence was the title, which a table of contents cannot see, so each is split
into a Heading 1 paragraph and a body paragraph.

The Contents is a Word TOC field rather than typed text, so page numbers are
correct after Word updates the field (it prompts on open, or Ctrl+A then F9).

Usage
-----
    python ops/docx/build_si_frontmatter.py FILE.docx --config sections.json
"""

from __future__ import annotations

import argparse
import json
import os
import re
import zipfile
from pathlib import Path
from xml.sax.saxutils import escape

DOC = "word/document.xml"


def para(text: str, *, style: str | None = None, bold: bool = False,
         size: int | None = None, align: str | None = None,
         color: str | None = None, italic: bool = False) -> str:
    ppr = "<w:pPr>"
    if style:
        ppr += f'<w:pStyle w:val="{style}"/>'
    if align:
        ppr += f'<w:jc w:val="{align}"/>'
    ppr += "</w:pPr>"
    rpr = "<w:rPr>"
    if bold:
        rpr += "<w:b/>"
    if italic:
        rpr += "<w:i/>"
    if color:
        rpr += f'<w:color w:val="{color}"/>'
    if size:
        rpr += f'<w:sz w:val="{size}"/><w:szCs w:val="{size}"/>'
    rpr += "</w:rPr>"
    return (f"<w:p>{ppr}<w:r>{rpr}"
            f'<w:t xml:space="preserve">{escape(text)}</w:t></w:r></w:p>')


def page_break() -> str:
    return '<w:p><w:r><w:br w:type="page"/></w:r></w:p>'


def toc_field() -> str:
    return (
        '<w:p><w:r><w:fldChar w:fldCharType="begin" w:dirty="true"/></w:r>'
        '<w:r><w:instrText xml:space="preserve"> TOC \\o "1-1" \\h \\z \\u </w:instrText></w:r>'
        '<w:r><w:fldChar w:fldCharType="separate"/></w:r>'
        '<w:r><w:rPr><w:i/></w:rPr><w:t xml:space="preserve">'
        'Right-click and choose Update Field to generate the contents.'
        '</w:t></w:r>'
        '<w:r><w:fldChar w:fldCharType="end"/></w:r></w:p>')


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("docx", type=Path)
    ap.add_argument("--config", type=Path, required=True)
    args = ap.parse_args(argv)
    cfg = json.loads(args.config.read_text())

    zin = zipfile.ZipFile(args.docx)
    items = {n: zin.read(n) for n in zin.namelist()}
    infos = zin.infolist()
    zin.close()
    doc = items[DOC].decode("utf-8")

    # split "Title. Body text..." into a Heading 1 and a body paragraph
    for sec in cfg["sections"]:
        anchor = sec["anchor"]
        m = re.search(r"<w:p\b[^>]*>(?:(?!</w:p>).)*?" + re.escape(anchor)
                      + r".*?</w:p>", doc, re.S)
        if not m:
            raise SystemExit(f"section paragraph not found: {anchor!r}")
        block = m.group(0)
        body = re.sub(r"^" + re.escape(anchor) + r"\s*", "", sec["full_text"])
        new = (para(sec["title"], style="Heading1", color=sec.get("color"))
               + re.sub(r'<w:t xml:space="preserve">.*?</w:t>',
                        f'<w:t xml:space="preserve">{escape(body)}</w:t>',
                        block, count=1, flags=re.S))
        doc = doc[:m.start()] + new + doc[m.end():]

    # figures heading becomes a section so it appears in the contents
    fh = cfg.get("figures_heading")
    if fh:
        m = re.search(r"<w:p\b[^>]*>(?:(?!</w:p>).)*?" + re.escape(fh["anchor"])
                      + r".*?</w:p>", doc, re.S)
        if m:
            doc = (doc[:m.start()]
                   + para(fh["title"], style="Heading1", color=fh.get("color"))
                   + doc[m.end():])

    # front matter, inserted at the very start of the body
    t = cfg["title_page"]
    front = [
        para("Supporting Information", style="Title", align="center"),
        para(t["manuscript_title"], bold=True, align="center", size=28),
        para("", align="center"),
        para(t["authors"], align="center", size=24),
        para(t["affiliation"], italic=True, align="center", size=22),
        para(t["correspondence"], align="center", size=22),
        page_break(),
        para("Contents", style="Heading1"),
        toc_field(),
        page_break(),
    ]
    bm = re.search(r"<w:body>", doc)
    doc = doc[:bm.end()] + "".join(front) + doc[bm.end():]

    items[DOC] = doc.encode("utf-8")
    tmp = args.docx.with_suffix(".tmp")
    with zipfile.ZipFile(tmp, "w", zipfile.ZIP_DEFLATED) as z:
        for info in infos:
            z.writestr(info, items[info.filename])
    os.replace(tmp, args.docx)
    print(f"{args.docx.name}: title page, contents and "
          f"{len(cfg['sections'])} section headings")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
