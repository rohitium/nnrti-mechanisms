#!/usr/bin/env python3
"""Edit a .docx in place: replace or append paragraph text without rebuilding.

The revision documents are edited by hand in Word between automated passes, so
regenerating them from a script would destroy that work. This unpacks the file
as it currently stands, changes only the paragraphs asked for, and repacks.

Operations
----------
replace-paragraph   Replace the runs of the paragraph containing MATCH.
append-to-paragraph Append a run to the paragraph containing MATCH.
insert-after        Add a NEW paragraph directly after the one containing MATCH.

Text added by this tool is red by default, matching the revision convention
that every automated addition is visible as such.

Usage
-----
    python ops/docx/edit_in_place.py replace-paragraph FILE.docx \
        --match "[PENDING - autocorrelation" --text "New text here."
    python ops/docx/edit_in_place.py append-to-paragraph FILE.docx \
        --match "some anchor" --text " Added sentence." --color 000000
"""

from __future__ import annotations

import argparse
import re
import shutil
import subprocess
import sys
import tempfile
import zipfile
from pathlib import Path
from xml.sax.saxutils import escape

DOC = "word/document.xml"


def _runs(text: str, color: str | None, bold: bool) -> str:
    rpr = ""
    if bold or color:
        rpr = "<w:rPr>" + ("<w:b/>" if bold else "") + \
              (f'<w:color w:val="{color}"/>' if color else "") + "</w:rPr>"
    return f'<w:r>{rpr}<w:t xml:space="preserve">{escape(text)}</w:t></w:r>'


def _para_text(p: str) -> str:
    return "".join(re.findall(r"<w:t[^>]*>(.*?)</w:t>", p, re.S))


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("op", choices=("replace-paragraph", "append-to-paragraph",
                                   "insert-after"))
    ap.add_argument("docx", type=Path)
    ap.add_argument("--match", required=True, help="substring identifying the paragraph")
    ap.add_argument("--text", required=True)
    ap.add_argument("--color", default="FF0000", help="hex, or 'none' for inherited")
    ap.add_argument("--bold", action="store_true")
    ap.add_argument("--occurrence", type=int, default=0,
                    help="which matching paragraph (0-based) when several match")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args(argv)

    color = None if args.color.lower() == "none" else args.color
    work = Path(tempfile.mkdtemp(prefix="docxedit_"))
    try:
        with zipfile.ZipFile(args.docx) as z:
            z.extractall(work)
        x = (work / DOC).read_text(encoding="utf-8")

        # <w:p .../> is an empty paragraph and is complete on its own. Without
        # the first alternative, a non-greedy match starting at one runs on to
        # the next </w:p> and returns two paragraphs glued together; editing
        # that emits an orphan run and an unmatched </w:p>.
        PARA_RE = r"<w:p\b[^>]*/>|<w:p\b[^>]*>.*?</w:p>"
        paras = re.findall(PARA_RE, x, re.S)
        hits = [p for p in paras if args.match in _para_text(p)]
        if not hits:
            print(f"no paragraph contains {args.match!r}", file=sys.stderr)
            return 2
        if args.occurrence >= len(hits):
            print(f"{len(hits)} matches, occurrence {args.occurrence} out of range", file=sys.stderr)
            return 2
        target = hits[args.occurrence]
        if x.count(target) != 1:
            print("target paragraph is not unique in the document", file=sys.stderr)
            return 2

        # Word writes an empty paragraph as a self-closing <w:p .../>. Reusing
        # that as an opening tag emits an orphan run and an unmatched </w:p>,
        # which is a mismatched-tag error Word will not open.
        if target.endswith("/>"):
            open_tag = "<w:p>"          # empty source paragraph; start a fresh one
        else:
            ppr = re.match(r"<w:p\b[^>]*>(<w:pPr>.*?</w:pPr>)?", target, re.S)
            open_tag = target[:ppr.end()]
        if args.op == "replace-paragraph":
            new = open_tag + _runs(args.text, color, args.bold) + "</w:p>"
        elif args.op == "insert-after":
            # A fresh paragraph carrying the matched paragraph's own properties,
            # so an inserted note inherits the surrounding body style rather
            # than Word's default.
            new = target + open_tag + _runs(args.text, color, args.bold) + "</w:p>"
        else:
            new = target[: -len("</w:p>")] + _runs(args.text, color, args.bold) + "</w:p>"

        print(f"matched {len(hits)} paragraph(s); editing #{args.occurrence}")
        print(f"  before: {_para_text(target)[:110]}")
        print(f"  after : {_para_text(new)[:110]}")
        if args.dry_run:
            return 0

        (work / DOC).write_text(x.replace(target, new), encoding="utf-8")
        tmp = args.docx.with_suffix(".tmp.docx")
        if tmp.exists():
            tmp.unlink()
        subprocess.run(["zip", "-Xrq", str(tmp.resolve()), "."], cwd=work, check=True)
        shutil.move(str(tmp), str(args.docx))
        print(f"  wrote {args.docx}")
        return 0
    finally:
        shutil.rmtree(work, ignore_errors=True)


if __name__ == "__main__":
    raise SystemExit(main())
