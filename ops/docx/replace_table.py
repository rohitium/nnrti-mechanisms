#!/usr/bin/env python3
"""Replace an existing table in a .docx with one built from a CSV.

insert_table.py only adds. Replacing a table that is already in the document
needs the old <w:tbl> removed first, which this does by matching text in its
header row rather than by position, so it cannot hit the wrong table.

Usage
-----
    python ops/docx/replace_table.py FILE.docx --header-contains "Single-point" \
        --csv new.csv [--color FF0000] [--font-size 14] [--widths ...]
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import sys
import zipfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from insert_table import build_table          # noqa: E402

DOC = "word/document.xml"


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("docx", type=Path)
    ap.add_argument("--header-contains", required=True,
                    help="text unique to the header row of the table to replace")
    ap.add_argument("--csv", type=Path, required=True)
    ap.add_argument("--color", default=None)
    ap.add_argument("--font-size", type=int, default=14)
    ap.add_argument("--widths", default=None)
    args = ap.parse_args(argv)

    with args.csv.open() as fh:
        rows = [r for r in csv.reader(fh) if any(c.strip() for c in r)]
    if not rows:
        raise SystemExit(f"{args.csv} has no rows")

    zin = zipfile.ZipFile(args.docx)
    items = {n: zin.read(n) for n in zin.namelist()}
    infos = zin.infolist()
    zin.close()
    doc = items[DOC].decode("utf-8")

    hits = []
    for m in re.finditer(r"<w:tbl>.*?</w:tbl>", doc, re.S):
        head = "".join(re.findall(r"<w:t(?: [^>]*)?>(.*?)</w:t>",
                                  m.group(0)[:4000], re.S))
        if args.header_contains in head:
            hits.append(m)
    if len(hits) != 1:
        raise SystemExit(f"{len(hits)} tables match {args.header_contains!r}; "
                         "need exactly one")
    old = hits[0]

    widths = [int(w) for w in args.widths.split(",")] if args.widths else None
    new = build_table(rows, color=args.color, sz=args.font_size, widths=widths)
    # build_table appends a trailing empty paragraph; the old table already has
    # one after it, so drop the duplicate
    new = new[:new.rindex("<w:p>")] if new.rstrip().endswith("</w:p>") else new

    doc = doc[:old.start()] + new + doc[old.end():]
    items[DOC] = doc.encode("utf-8")

    tmp = args.docx.with_suffix(".tmp")
    with zipfile.ZipFile(tmp, "w", zipfile.ZIP_DEFLATED) as z:
        for info in infos:
            z.writestr(info, items[info.filename])
    os.replace(tmp, args.docx)
    print(f"{args.docx.name}: replaced table ({len(rows)}x{max(len(r) for r in rows)}), "
          f"tables now {doc.count('<w:tbl>')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
