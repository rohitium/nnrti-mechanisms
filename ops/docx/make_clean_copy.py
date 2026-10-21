#!/usr/bin/env python3
"""Produce the clean, unmarked copy from the red-marked revision.

ACS asks for two manuscript files: a clean one for publication and a marked one
showing the changes, uploaded as Supporting Information for Review. Our
revisions are marked by colouring every addition red, so the clean copy is the
same document with that colouring removed.

This strips only the red (FF0000) colour runs, leaving the text, and leaves all
other formatting untouched. It does not remove content: the clean copy and the
marked copy say the same thing, which is the point -- a clean copy that
differed in substance would be a different paper.

Usage
-----
    python ops/docx/make_clean_copy.py MARKED.docx CLEAN.docx
"""

from __future__ import annotations

import re
import shutil
import subprocess
import sys
import tempfile
import zipfile
from pathlib import Path

DOC = "word/document.xml"
RED = re.compile(r'<w:color w:val="FF0000"\s*/>', re.I)


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if len(argv) != 2:
        print(__doc__)
        return 2
    marked, clean = Path(argv[0]), Path(argv[1])

    work = Path(tempfile.mkdtemp(prefix="cleancopy_"))
    try:
        with zipfile.ZipFile(marked) as z:
            z.extractall(work)
        x = (work / DOC).read_text(encoding="utf-8")
        n = len(RED.findall(x))
        x = RED.sub("", x)
        # An <w:rPr> left with no children is legal but pointless; drop the empties.
        x = re.sub(r"<w:rPr>\s*</w:rPr>", "", x)
        (work / DOC).write_text(x, encoding="utf-8")

        if clean.exists():
            clean.unlink()
        subprocess.run(["zip", "-Xrq", str(clean.resolve()), "."], cwd=work, check=True)
        print(f"removed {n} red colour tags -> {clean}")
        return 0
    finally:
        shutil.rmtree(work, ignore_errors=True)


if __name__ == "__main__":
    raise SystemExit(main())
