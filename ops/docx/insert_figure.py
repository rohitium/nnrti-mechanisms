#!/usr/bin/env python3
"""Put an image above a figure caption that currently has none.

Supplementary Figures 3-5 existed as captions with nothing behind them. This
wires an image into the package -- media file, relationship, content type --
and inserts the drawing in its own paragraph immediately before the caption.

Width is given in inches and the height follows the image's aspect ratio, so a
figure cannot be distorted by a wrong number.

Usage
-----
    python ops/docx/insert_figure.py FILE.docx IMAGE.png \
        --before "Supplementary Figure 3." [--width-in 6.5]
"""

from __future__ import annotations

import argparse
import os
import re
import zipfile
from pathlib import Path

DOC = "word/document.xml"
RELS = "word/_rels/document.xml.rels"
CT = "[Content_Types].xml"
EMU = 914400


def main(argv: list[str] | None = None) -> int:
    from PIL import Image

    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("docx", type=Path)
    ap.add_argument("image", type=Path)
    ap.add_argument("--before", required=True, help="text of the caption paragraph")
    ap.add_argument("--width-in", type=float, default=6.5)
    ap.add_argument("--name", default=None)
    args = ap.parse_args(argv)

    zin = zipfile.ZipFile(args.docx)
    items = {n: zin.read(n) for n in zin.namelist()}
    infos = zin.infolist()
    zin.close()
    doc = items[DOC].decode("utf-8")

    name = args.name or args.image.stem
    if f'name="{name}"' in doc:
        raise SystemExit(f"a drawing named {name!r} is already in the document")

    ext = args.image.suffix.lstrip(".").lower()
    media = f"word/media/{name}.{ext}"
    items[media] = args.image.read_bytes()

    rels = items[RELS].decode("utf-8")
    rid = max(int(m) for m in re.findall(r'Id="rId(\d+)"', rels)) + 1
    rels = rels.replace(
        "</Relationships>",
        f'<Relationship Id="rId{rid}" Type="http://schemas.openxmlformats.org/'
        f'officeDocument/2006/relationships/image" Target="media/{name}.{ext}"/>'
        "</Relationships>")
    items[RELS] = rels.encode("utf-8")

    ct = items[CT].decode("utf-8")
    if f'Extension="{ext}"' not in ct:
        mime = {"png": "image/png", "jpg": "image/jpeg",
                "jpeg": "image/jpeg", "tif": "image/tiff"}.get(ext, "image/png")
        ct = ct.replace("<Default", f'<Default Extension="{ext}" '
                        f'ContentType="{mime}"/><Default', 1)
        items[CT] = ct.encode("utf-8")

    with Image.open(args.image) as im:
        w, h = im.size
    cx = int(args.width_in * EMU)
    cy = int(cx * h / w)

    drawing = (
        f'<w:p><w:pPr><w:jc w:val="center"/></w:pPr><w:r><w:drawing>'
        f'<wp:inline distT="0" distB="0" distL="0" distR="0" '
        f'xmlns:wp="http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing">'
        f'<wp:extent cx="{cx}" cy="{cy}"/><wp:docPr id="{900 + rid}" name="{name}"/>'
        f'<a:graphic xmlns:a="http://schemas.openxmlformats.org/drawingml/2006/main">'
        f'<a:graphicData uri="http://schemas.openxmlformats.org/drawingml/2006/picture">'
        f'<pic:pic xmlns:pic="http://schemas.openxmlformats.org/drawingml/2006/picture">'
        f'<pic:nvPicPr><pic:cNvPr id="{900 + rid}" name="{name}"/><pic:cNvPicPr/></pic:nvPicPr>'
        f'<pic:blipFill><a:blip xmlns:r="http://schemas.openxmlformats.org/'
        f'officeDocument/2006/relationships" r:embed="rId{rid}"/>'
        f'<a:stretch><a:fillRect/></a:stretch></pic:blipFill>'
        f'<pic:spPr><a:xfrm><a:off x="0" y="0"/><a:ext cx="{cx}" cy="{cy}"/></a:xfrm>'
        f'<a:prstGeom prst="rect"><a:avLst/></a:prstGeom></pic:spPr>'
        f'</pic:pic></a:graphicData></a:graphic></wp:inline></w:drawing></w:r></w:p>')

    pat = re.compile(r"<w:p\b[^>]*>(?:(?!</w:p>).)*?" + re.escape(args.before)
                     + r".*?</w:p>", re.S)
    m = pat.search(doc)
    if not m:
        raise SystemExit(f"caption not found: {args.before!r}")
    doc = doc[:m.start()] + drawing + doc[m.start():]
    items[DOC] = doc.encode("utf-8")

    tmp = args.docx.with_suffix(".tmp")
    with zipfile.ZipFile(tmp, "w", zipfile.ZIP_DEFLATED) as z:
        for info in infos:
            z.writestr(info, items[info.filename])
        for extra in (media,):
            if extra not in {i.filename for i in infos}:
                z.writestr(extra, items[extra])
    os.replace(tmp, args.docx)
    print(f"inserted {name} at {args.width_in}x{cy/EMU:.2f} in before "
          f"{args.before!r}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
