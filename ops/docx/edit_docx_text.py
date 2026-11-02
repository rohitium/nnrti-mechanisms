#!/usr/bin/env python3
"""Replace exact strings inside a .docx without disturbing anything else.

Word splits a sentence across many runs, often mid-word, so a naive search of
the XML misses most real text and a naive rewrite destroys formatting. This
walks each paragraph, reconstructs its visible text, finds the target there,
and edits only the runs the match actually touches. Every other run, and every
field, bookmark, comment and style, is left untouched.

Guard rails, because this edits a document nobody wants to retype:

  * a replacement that does not match exactly once in the file is refused, so
    an ambiguous or already-applied edit cannot be half-made;
  * paragraphs containing a field (Zotero citations) are skipped by default,
    since editing around a field can orphan it;
  * the formatting of the run where the match starts is what the replacement
    inherits, so text in a red "changed" run stays red.

Usage
-----
    python ops/docx/edit_docx_text.py IN.docx OUT.docx edits.json
"""

from __future__ import annotations

import json
import shutil
import sys
import zipfile
from pathlib import Path
import xml.etree.ElementTree as ET

W = "{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"
ET.register_namespace("w", W[1:-1])


def _runs(par):
    return [r for r in par.iter(W + "r")]


def _run_text(run) -> str:
    return "".join(t.text or "" for t in run.iter(W + "t"))


def _set_run_text(run, text: str) -> None:
    ts = list(run.iter(W + "t"))
    if not ts:
        return
    ts[0].text = text
    # xml:space so leading/trailing spaces survive the round trip.
    ts[0].set("{http://www.w3.org/XML/1998/namespace}space", "preserve")
    for extra in ts[1:]:
        extra.text = ""


def replace_in_paragraph(par, old: str, new: str) -> bool:
    runs = _runs(par)
    texts = [_run_text(r) for r in runs]
    joined = "".join(texts)
    idx = joined.find(old)
    if idx < 0:
        return False

    # Map the match onto the runs it spans.
    start, end = idx, idx + len(old)
    pos = 0
    first = None
    for i, t in enumerate(texts):
        nxt = pos + len(t)
        if first is None and nxt > start:
            first = i
            head = t[: start - pos]
        if pos < end <= nxt:
            tail = t[end - pos :]
            last = i
            break
        pos = nxt
    else:
        return False

    if first == last:
        _set_run_text(runs[first], head + new + tail)
    else:
        _set_run_text(runs[first], head + new)
        for i in range(first + 1, last):
            _set_run_text(runs[i], "")
        _set_run_text(runs[last], tail)
    return True


def has_field(par) -> bool:
    x = ET.tostring(par, encoding="unicode")
    return "fldChar" in x or "ZOTERO_ITEM" in x or "instrText" in x


def apply_edits(src: Path, dest: Path, edits: list[dict],
                allow_fields: bool = False) -> None:
    shutil.copyfile(src, dest)
    with zipfile.ZipFile(src) as z:
        names = z.namelist()
        blobs = {n: z.read(n) for n in names}

    root = ET.fromstring(blobs["word/document.xml"])
    body = root.find(f"{W}body")
    pars = [p for p in body.iter(f"{W}p")]

    report = []
    for e in edits:
        old, new = e["old"], e["new"]
        hits = [p for p in pars
                if old in "".join(_run_text(r) for r in _runs(p))
                and (allow_fields or not has_field(p))]
        if len(hits) != 1:
            raise SystemExit(
                f"REFUSED: {len(hits)} matches for {old[:70]!r}; "
                "an edit must match exactly one paragraph")
        if not replace_in_paragraph(hits[0], old, new):
            raise SystemExit(f"REFUSED: could not place {old[:70]!r}")
        report.append((old, new))

    blobs["word/document.xml"] = ET.tostring(root, encoding="UTF-8",
                                             xml_declaration=True)
    with zipfile.ZipFile(dest, "w", zipfile.ZIP_DEFLATED) as z:
        for n in names:
            z.writestr(n, blobs[n])

    print(f"{len(report)} edits applied -> {dest}")
    for old, new in report:
        print(f"   - {old[:68]}")
        print(f"   + {new[:68]}")


def main() -> int:
    if len(sys.argv) != 4:
        print(__doc__)
        return 2
    src, dest, spec = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
    apply_edits(src, dest, json.loads(spec.read_text()))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
