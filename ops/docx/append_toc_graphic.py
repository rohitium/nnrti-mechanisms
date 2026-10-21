#!/usr/bin/env python3
"""Append the TOC graphic to the last page of a manuscript .docx.

ACS asks for the TOC graphic on the last page, labelled "TOC Graphic", with no
caption. This adds a page break, the label, and the image, wiring up the media
file, the relationship and the content-type override that Word needs.

Idempotent: if a "TOC Graphic" label is already present the file is left alone,
so re-running after other edits does not stack duplicates.

Usage
-----
    python ops/docx/append_toc_graphic.py MANUSCRIPT.docx IMAGE.png [--width-in 3.3]
"""

from __future__ import annotations

import argparse
import re
import shutil
import subprocess
import tempfile
import zipfile
from pathlib import Path

EMU_PER_INCH = 914400
DOC = "word/document.xml"
RELS = "word/_rels/document.xml.rels"
CT = "[Content_Types].xml"


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("docx", type=Path)
    ap.add_argument("image", type=Path)
    ap.add_argument("--width-in", type=float, default=3.3)
    ap.add_argument("--height-in", type=float, default=1.8)
    args = ap.parse_args(argv)

    work = Path(tempfile.mkdtemp(prefix="toc_"))
    try:
        with zipfile.ZipFile(args.docx) as z:
            z.extractall(work)
        doc = (work / DOC).read_text(encoding="utf-8")
        if "TOC Graphic" in doc:
            print("already present; nothing to do")
            return 0

        ext = args.image.suffix.lstrip(".").lower()
        media = work / "word" / "media"
        media.mkdir(parents=True, exist_ok=True)
        shutil.copy(args.image, media / f"toc_graphic.{ext}")

        rels = (work / RELS).read_text(encoding="utf-8")
        used = {int(m) for m in re.findall(r'Id="rId(\d+)"', rels)}
        rid = f"rId{max(used) + 1 if used else 1}"
        rels = rels.replace("</Relationships>",
            f'<Relationship Id="{rid}" Type="http://schemas.openxmlformats.org/'
            f'officeDocument/2006/relationships/image" Target="media/toc_graphic.{ext}"/>'
            "</Relationships>")
        (work / RELS).write_text(rels, encoding="utf-8")

        ct = (work / CT).read_text(encoding="utf-8")
        if f'Extension="{ext}"' not in ct:
            mime = {"png": "image/png", "tif": "image/tiff", "tiff": "image/tiff",
                    "jpg": "image/jpeg", "jpeg": "image/jpeg"}.get(ext, "image/png")
            ct = ct.replace("</Types>",
                            f'<Default Extension="{ext}" ContentType="{mime}"/></Types>')
            (work / CT).write_text(ct, encoding="utf-8")

        cx, cy = int(args.width_in * EMU_PER_INCH), int(args.height_in * EMU_PER_INCH)
        block = (
            '<w:p><w:r><w:br w:type="page"/></w:r></w:p>'
            '<w:p><w:pPr><w:spacing w:after="240"/></w:pPr>'
            '<w:r><w:rPr><w:b/></w:rPr><w:t>TOC Graphic</w:t></w:r></w:p>'
            '<w:p><w:r><w:drawing>'
            f'<wp:inline distT="0" distB="0" distL="0" distR="0">'
            f'<wp:extent cx="{cx}" cy="{cy}"/><wp:docPr id="1001" name="TOC Graphic"/>'
            '<a:graphic xmlns:a="http://schemas.openxmlformats.org/drawingml/2006/main">'
            '<a:graphicData uri="http://schemas.openxmlformats.org/drawingml/2006/picture">'
            '<pic:pic xmlns:pic="http://schemas.openxmlformats.org/drawingml/2006/picture">'
            '<pic:nvPicPr><pic:cNvPr id="1001" name="TOC Graphic"/><pic:cNvPicPr/></pic:nvPicPr>'
            f'<pic:blipFill><a:blip r:embed="{rid}"/><a:stretch><a:fillRect/></a:stretch></pic:blipFill>'
            '<pic:spPr><a:xfrm><a:off x="0" y="0"/>'
            f'<a:ext cx="{cx}" cy="{cy}"/></a:xfrm>'
            '<a:prstGeom prst="rect"><a:avLst/></a:prstGeom></pic:spPr>'
            '</pic:pic></a:graphicData></a:graphic></wp:inline>'
            '</w:drawing></w:r></w:p>'
        )
        # Insert before the BODY-level sectPr, which is the last one in the file.
        # A lazy regex from the first <w:sectPr> matches a section break inside a
        # paragraph's properties instead, and placing a <w:p> there is invalid.
        body_end = doc.rindex("</w:body>")
        sect = doc.rfind("<w:sectPr", 0, body_end)
        at = sect if sect != -1 else body_end
        doc = doc[:at] + block + doc[at:]

        # The drawing markup needs the DrawingML namespaces; Word documents
        # without an inline image do not declare them.
        m = re.search(r"<w:document\b[^>]*>", doc)
        root = m.group(0)
        add = ""
        if "xmlns:a=" not in root:
            add += ' xmlns:a="http://schemas.openxmlformats.org/drawingml/2006/main"'
        if "xmlns:pic=" not in root:
            add += ' xmlns:pic="http://schemas.openxmlformats.org/drawingml/2006/picture"'
        if add:
            doc = doc[:m.end() - 1] + add + doc[m.end() - 1:]
        (work / DOC).write_text(doc, encoding="utf-8")

        tmp = args.docx.with_suffix(".tmp.docx")
        if tmp.exists():
            tmp.unlink()
        subprocess.run(["zip", "-Xrq", str(tmp.resolve()), "."], cwd=work, check=True)
        shutil.move(str(tmp), str(args.docx))
        print(f"appended TOC graphic ({args.width_in}x{args.height_in} in) to {args.docx}")
        return 0
    finally:
        shutil.rmtree(work, ignore_errors=True)


if __name__ == "__main__":
    raise SystemExit(main())
