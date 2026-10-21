#!/usr/bin/env python3
"""Insert a Word table into an existing .docx, in place.

The revision documents are edited by hand in Word between automated passes,
so they are never rebuilt. This unpacks the file as it stands, splices one
<w:tbl> in after the paragraph containing MATCH, and repacks.

Tables are written with a plain single-line border on every cell, a bold
header row, and a caption paragraph is NOT added -- captions already exist in
these documents as their own paragraphs, and the table goes after them.

Usage
-----
    python ops/docx/insert_table.py FILE.docx --match "Supplementary Table 5." \
        --csv table.csv [--color FF0000] [--font-size 16] [--widths 2000,1200,...]

--color applies to the body text (header stays black); use it to mark a table
as an addition under the revision convention. --font-size is in half-points,
so 16 means 8 pt.
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import zipfile
from pathlib import Path
from xml.sax.saxutils import escape

DOC = "word/document.xml"


def _cell(text: str, *, bold: bool, color: str | None, sz: int, width: int) -> str:
    rpr = "<w:rPr>"
    if bold:
        rpr += "<w:b/>"
    if color:
        rpr += f'<w:color w:val="{color}"/>'
    rpr += f'<w:sz w:val="{sz}"/><w:szCs w:val="{sz}"/></w:rPr>'
    # keep line breaks inside a cell as separate runs
    runs = f'<w:r>{rpr}<w:t xml:space="preserve">{escape(text)}</w:t></w:r>'
    return (
        f'<w:tc><w:tcPr><w:tcW w:w="{width}" w:type="dxa"/>'
        f'<w:vAlign w:val="center"/></w:tcPr>'
        f'<w:p><w:pPr><w:spacing w:after="0" w:line="240" w:lineRule="auto"/>'
        f'<w:jc w:val="center"/></w:pPr>{runs}</w:p></w:tc>'
    )


def build_table(rows: list[list[str]], *, color: str | None, sz: int,
                widths: list[int] | None) -> str:
    ncol = max(len(r) for r in rows)
    if widths and len(widths) != ncol:
        raise SystemExit(f"--widths has {len(widths)} entries, table has {ncol} columns")
    widths = widths or [int(9360 / ncol)] * ncol
    borders = "".join(
        f'<w:{e} w:val="single" w:sz="4" w:space="0" w:color="auto"/>'
        for e in ("top", "left", "bottom", "right", "insideH", "insideV"))
    xml = ['<w:tbl><w:tblPr><w:tblStyle w:val="TableGrid"/>'
           '<w:tblW w:w="0" w:type="auto"/>'
           f'<w:tblBorders>{borders}</w:tblBorders></w:tblPr>']
    xml.append("<w:tblGrid>" + "".join(f'<w:gridCol w:w="{w}"/>' for w in widths)
               + "</w:tblGrid>")
    for i, row in enumerate(rows):
        head = i == 0
        cells = "".join(
            _cell(row[j] if j < len(row) else "", bold=head,
                  color=None if head else color, sz=sz, width=widths[j])
            for j in range(ncol))
        trpr = '<w:trPr><w:tblHeader/></w:trPr>' if head else ""
        xml.append(f"<w:tr>{trpr}{cells}</w:tr>")
    xml.append("</w:tbl>")
    # Word requires a paragraph after a table; without it two adjacent tables
    # merge and the document fails to open.
    xml.append('<w:p><w:pPr><w:spacing w:after="0"/></w:pPr></w:p>')
    return "".join(xml)


def insert_after(doc: str, match: str, tbl: str, occurrence: int = 0) -> str:
    paras = [m for m in re.finditer(r"<w:p[ >].*?</w:p>", doc, re.S)]
    hits = [m for m in paras
            if match in html_text(m.group(0))]
    if not hits:
        raise SystemExit(f"no paragraph contains {match!r}")
    if occurrence >= len(hits):
        raise SystemExit(f"{len(hits)} paragraphs match {match!r}; "
                         f"--occurrence {occurrence} is out of range")
    at = hits[occurrence].end()
    return doc[:at] + tbl + doc[at:]


def html_text(p: str) -> str:
    import html as _h
    p = re.sub(r"<w:instrText[^>]*>.*?</w:instrText>", "", p, flags=re.S)
    return _h.unescape("".join(re.findall(r"<w:t[^>]*>(.*?)</w:t>", p, re.S)))


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("docx", type=Path)
    ap.add_argument("--match", required=True)
    ap.add_argument("--csv", type=Path, required=True)
    ap.add_argument("--color", default=None, help="hex for body text, e.g. FF0000")
    ap.add_argument("--font-size", type=int, default=16, help="half-points (16 = 8pt)")
    ap.add_argument("--widths", default=None, help="comma-separated dxa widths")
    ap.add_argument("--occurrence", type=int, default=0)
    args = ap.parse_args(argv)

    with args.csv.open() as fh:
        rows = [r for r in csv.reader(fh) if any(c.strip() for c in r)]
    if not rows:
        raise SystemExit(f"{args.csv} has no rows")

    widths = [int(w) for w in args.widths.split(",")] if args.widths else None
    tbl = build_table(rows, color=args.color, sz=args.font_size, widths=widths)

    zin = zipfile.ZipFile(args.docx)
    items = {n: zin.read(n) for n in zin.namelist()}
    infos = zin.infolist()
    zin.close()
    doc = items[DOC].decode("utf8")
    before = doc.count("<w:tbl>")
    doc = insert_after(doc, args.match, tbl, args.occurrence)
    items[DOC] = doc.encode("utf8")

    tmp = args.docx.with_suffix(".tmp")
    with zipfile.ZipFile(tmp, "w", zipfile.ZIP_DEFLATED) as z:
        for info in infos:
            z.writestr(info, items[info.filename])
    os.replace(tmp, args.docx)
    print(f"{args.docx.name}: inserted {len(rows)}x{max(len(r) for r in rows)} table "
          f"after {args.match!r} (tables {before} -> {doc.count('<w:tbl>')})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
