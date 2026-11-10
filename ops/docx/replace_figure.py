#!/usr/bin/env python3
"""Swap the image that already sits above a figure caption.

insert_figure.py wires a new image into a document that has none. This is the
other half: the caption already has a drawing, and the art behind it needs to
change. Replacing the bytes of word/media/imageN.png alone is not enough --
the drawing carries its own wp:extent and a:ext in EMU, so a new image with a
different aspect ratio renders stretched with no warning anywhere.

The width is preserved from the existing drawing unless --width-in is given,
and the height is recomputed from the new image's aspect ratio, so a figure
cannot be distorted by a stale extent.

Usage
-----
    python ops/docx/replace_figure.py FILE.docx IMAGE.png \
        --before "Supplementary Figure 6." [--width-in 6.5] [--out OUT.docx]
"""

from __future__ import annotations

import argparse
import re
import shutil
import zipfile
from pathlib import Path

DOC = "word/document.xml"
RELS = "word/_rels/document.xml.rels"
EMU = 914400


def main(argv: list[str] | None = None) -> int:
    from PIL import Image

    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("docx", type=Path)
    ap.add_argument("image", type=Path)
    ap.add_argument("--before", required=True,
                    help="leading text of the caption paragraph the image sits above")
    ap.add_argument("--width-in", type=float, default=None)
    ap.add_argument("--out", type=Path, default=None)
    args = ap.parse_args(argv)

    zin = zipfile.ZipFile(args.docx)
    items = {n: zin.read(n) for n in zin.namelist()}
    infos = zin.infolist()
    zin.close()

    doc = items[DOC].decode("utf8")
    rels = items[RELS].decode("utf8")
    rid_target = {m.group(1): m.group(2) for m in
                  re.finditer(r'Id="(rId\d+)"[^>]*Target="([^"]+)"', rels)}

    paras = list(re.finditer(r"<w:p[ >].*?</w:p>", doc, re.S))

    def visible(xml: str) -> str:
        return re.sub(r"<[^>]+>", "", "".join(
            re.findall(r"<w:t[^>]*>(.*?)</w:t>", xml, re.S)))

    cap = [i for i, m in enumerate(paras)
           if visible(m.group(0)).startswith(args.before)]
    if len(cap) != 1:
        raise SystemExit(f"caption {args.before!r} matched {len(cap)} paragraphs; "
                         "it must match exactly one")
    # The drawing may be in the caption's own paragraph -- which is how
    # Figure 5 is laid out in the manuscript, art and caption text in one
    # <w:p> -- or in the nearest preceding paragraph that holds one, which is
    # how the Supplementary Figures are laid out. Checking the caption
    # paragraph first means one tool handles both, and a blank paragraph
    # between art and caption does not break the match.
    img_par = None
    if "r:embed=" in paras[cap[0]].group(0):
        img_par = cap[0]
    else:
        for i in range(cap[0] - 1, -1, -1):
            if "r:embed=" in paras[i].group(0):
                img_par = i
                break
            if visible(paras[i].group(0)).strip():
                break
    if img_par is None:
        raise SystemExit(f"no drawing found above {args.before!r}; "
                         "use ops/docx/insert_figure.py instead")

    par_xml = paras[img_par].group(0)
    embeds = re.findall(r'r:embed="(rId\d+)"', par_xml)
    if len(embeds) != 1:
        raise SystemExit(f"expected one drawing above the caption, found {len(embeds)}")
    target = rid_target[embeds[0]]
    media = "word/" + target.lstrip("/") if not target.startswith("word/") else target
    if media not in items:
        raise SystemExit(f"{media} not in the package")

    # Word writes a space before the self-closing slash and ElementTree does
    # not, so both spellings have to be accepted or the tool works on exactly
    # one of the two kinds of document in this repository.
    old_extent = re.search(r'<wp:extent cx="(\d+)" cy="(\d+)"\s*/>', par_xml)
    if not old_extent:
        raise SystemExit("drawing has no wp:extent")
    cx_old = int(old_extent.group(1))
    width_emu = int(round(args.width_in * EMU)) if args.width_in else cx_old

    with Image.open(args.image) as im:
        w, h = im.size
    height_emu = int(round(width_emu * h / w))

    new_par = re.sub(r'<wp:extent cx="\d+" cy="\d+"\s*/>',
                     f'<wp:extent cx="{width_emu}" cy="{height_emu}" />', par_xml)
    # Only the <a:ext> that carries cx/cy; <a:extLst> and the uri-keyed
    # <a:ext uri="..."> inside it must not be touched.
    new_par = re.sub(r'(<a:ext )cx="\d+" cy="\d+"',
                     rf'\1cx="{width_emu}" cy="{height_emu}"', new_par)
    # A cropped or resized drawing can also carry a stale srcRect; a full
    # replacement image is never meant to be cropped.
    new_par = re.sub(r'<a:srcRect[^/]*/>', '<a:srcRect/>', new_par)

    doc = doc[:paras[img_par].start()] + new_par + doc[paras[img_par].end():]
    items[DOC] = doc.encode("utf8")
    items[media] = args.image.read_bytes()

    out = args.out or args.docx
    if out == args.docx:
        shutil.copy2(args.docx, args.docx.with_suffix(args.docx.suffix + ".bak"))
    with zipfile.ZipFile(out, "w", zipfile.ZIP_DEFLATED) as zout:
        for info in infos:
            zout.writestr(info, items[info.filename])
    print(f"replaced {media} ({len(args.image.read_bytes())} bytes) above "
          f"{args.before!r}")
    print(f"extent {cx_old} -> {width_emu} x {height_emu} EMU "
          f"({width_emu/EMU:.2f} x {height_emu/EMU:.2f} in)")
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
