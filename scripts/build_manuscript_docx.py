#!/usr/bin/env python3
from __future__ import annotations

import argparse
import datetime as _dt
import re
import zipfile
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable
import xml.etree.ElementTree as ET


NS = {
    "w": "http://schemas.openxmlformats.org/wordprocessingml/2006/main",
    "r": "http://schemas.openxmlformats.org/officeDocument/2006/relationships",
    "wp": "http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing",
    "a": "http://schemas.openxmlformats.org/drawingml/2006/main",
    "pic": "http://schemas.openxmlformats.org/drawingml/2006/picture",
    "rels": "http://schemas.openxmlformats.org/package/2006/relationships",
    "cp": "http://schemas.openxmlformats.org/package/2006/metadata/core-properties",
    "dc": "http://purl.org/dc/elements/1.1/",
    "dcterms": "http://purl.org/dc/terms/",
    "dcmitype": "http://purl.org/dc/dcmitype/",
    "xsi": "http://www.w3.org/2001/XMLSchema-instance",
    "ep": "http://schemas.openxmlformats.org/officeDocument/2006/extended-properties",
    "vt": "http://schemas.openxmlformats.org/officeDocument/2006/docPropsVTypes",
}

for _prefix, _uri in NS.items():
    ET.register_namespace(_prefix, _uri)


def _qn(prefix: str, tag: str) -> str:
    return f"{{{NS[prefix]}}}{tag}"


EMU_PER_INCH = 914400
EMU_PER_PX_96DPI = EMU_PER_INCH // 96  # 9525


def _png_dimensions(path: Path) -> tuple[int, int]:
    data = path.read_bytes()
    if len(data) < 24 or data[:8] != b"\x89PNG\r\n\x1a\n":
        raise ValueError(f"Not a PNG file: {path}")
    width = int.from_bytes(data[16:20], "big")
    height = int.from_bytes(data[20:24], "big")
    if width <= 0 or height <= 0:
        raise ValueError(f"Invalid PNG dimensions in {path}: {width}x{height}")
    return width, height


@dataclass(frozen=True)
class _ImagePart:
    rid: str
    filename: str
    data: bytes
    cx: int
    cy: int
    docpr_id: int


class _DocxWriter:
    def __init__(self, title: str) -> None:
        self._title = title
        self._body = ET.Element(_qn("w", "body"))
        self._images: list[_ImagePart] = []
        self._next_rid = 1
        self._next_docpr = 1

    def _add_p_pr(
        self,
        p: ET.Element,
        align: str | None,
        spacing_before_twips: int | None,
        spacing_after_twips: int | None,
        indent_left_twips: int | None,
    ) -> None:
        if not any(v is not None for v in (align, spacing_before_twips, spacing_after_twips, indent_left_twips)):
            return
        ppr = ET.SubElement(p, _qn("w", "pPr"))
        if align is not None:
            jc = ET.SubElement(ppr, _qn("w", "jc"))
            jc.set(_qn("w", "val"), align)
        if spacing_before_twips is not None or spacing_after_twips is not None:
            sp = ET.SubElement(ppr, _qn("w", "spacing"))
            if spacing_before_twips is not None:
                sp.set(_qn("w", "before"), str(int(spacing_before_twips)))
            if spacing_after_twips is not None:
                sp.set(_qn("w", "after"), str(int(spacing_after_twips)))
        if indent_left_twips is not None:
            ind = ET.SubElement(ppr, _qn("w", "ind"))
            ind.set(_qn("w", "left"), str(int(indent_left_twips)))

    def _add_run(
        self,
        p: ET.Element,
        text: str,
        bold: bool = False,
        italic: bool = False,
        size_half_points: int | None = None,
    ) -> None:
        r = ET.SubElement(p, _qn("w", "r"))
        if bold or italic or size_half_points is not None:
            rpr = ET.SubElement(r, _qn("w", "rPr"))
            if bold:
                ET.SubElement(rpr, _qn("w", "b"))
            if italic:
                ET.SubElement(rpr, _qn("w", "i"))
            if size_half_points is not None:
                sz = ET.SubElement(rpr, _qn("w", "sz"))
                sz.set(_qn("w", "val"), str(int(size_half_points)))
                szcs = ET.SubElement(rpr, _qn("w", "szCs"))
                szcs.set(_qn("w", "val"), str(int(size_half_points)))
        t = ET.SubElement(r, _qn("w", "t"))
        # Preserve leading/trailing whitespace if present.
        if text[:1].isspace() or text[-1:].isspace():
            t.set("{http://www.w3.org/XML/1998/namespace}space", "preserve")
        t.text = text

    def add_paragraph(
        self,
        text: str,
        *,
        align: str | None = None,
        bold: bool = False,
        italic: bool = False,
        size_half_points: int | None = None,
        spacing_before_twips: int | None = None,
        spacing_after_twips: int | None = None,
        indent_left_twips: int | None = None,
    ) -> None:
        p = ET.SubElement(self._body, _qn("w", "p"))
        self._add_p_pr(p, align, spacing_before_twips, spacing_after_twips, indent_left_twips)
        self._add_run(p, text, bold=bold, italic=italic, size_half_points=size_half_points)

    def add_image(self, path: Path, caption: str | None) -> None:
        if not path.exists():
            raise FileNotFoundError(f"Missing image: {path}")
        if path.suffix.lower() != ".png":
            raise ValueError(f"Only PNG images are supported (got {path.suffix}): {path}")

        width_px, height_px = _png_dimensions(path)
        # Keep within printable width for Letter with 1in margins: 6.5 inches.
        max_width_in = 6.5
        width_in = width_px / 96.0
        scale = min(1.0, max_width_in / max(width_in, 1e-9))
        cx = int(width_px * scale * EMU_PER_PX_96DPI)
        cy = int(height_px * scale * EMU_PER_PX_96DPI)

        rid = f"rId{self._next_rid}"
        self._next_rid += 1
        docpr_id = self._next_docpr
        self._next_docpr += 1

        image_name = f"image{len(self._images) + 1}.png"
        self._images.append(
            _ImagePart(
                rid=rid,
                filename=image_name,
                data=path.read_bytes(),
                cx=cx,
                cy=cy,
                docpr_id=docpr_id,
            )
        )

        # Image paragraph (centered).
        p = ET.SubElement(self._body, _qn("w", "p"))
        self._add_p_pr(p, align="center", spacing_before_twips=120, spacing_after_twips=120, indent_left_twips=None)
        r = ET.SubElement(p, _qn("w", "r"))
        drawing = ET.SubElement(r, _qn("w", "drawing"))

        inline = ET.SubElement(drawing, _qn("wp", "inline"))
        inline.set("distT", "0")
        inline.set("distB", "0")
        inline.set("distL", "0")
        inline.set("distR", "0")

        extent = ET.SubElement(inline, _qn("wp", "extent"))
        extent.set("cx", str(cx))
        extent.set("cy", str(cy))

        docpr = ET.SubElement(inline, _qn("wp", "docPr"))
        docpr.set("id", str(docpr_id))
        docpr.set("name", f"Picture {docpr_id}")

        graphic = ET.SubElement(inline, _qn("a", "graphic"))
        graphic_data = ET.SubElement(graphic, _qn("a", "graphicData"))
        graphic_data.set("uri", "http://schemas.openxmlformats.org/drawingml/2006/picture")

        pic = ET.SubElement(graphic_data, _qn("pic", "pic"))

        nv = ET.SubElement(pic, _qn("pic", "nvPicPr"))
        cnvpr = ET.SubElement(nv, _qn("pic", "cNvPr"))
        cnvpr.set("id", "0")
        cnvpr.set("name", image_name)
        ET.SubElement(nv, _qn("pic", "cNvPicPr"))

        blip_fill = ET.SubElement(pic, _qn("pic", "blipFill"))
        blip = ET.SubElement(blip_fill, _qn("a", "blip"))
        blip.set(_qn("r", "embed"), rid)
        stretch = ET.SubElement(blip_fill, _qn("a", "stretch"))
        ET.SubElement(stretch, _qn("a", "fillRect"))

        sppr = ET.SubElement(pic, _qn("pic", "spPr"))
        xfrm = ET.SubElement(sppr, _qn("a", "xfrm"))
        off = ET.SubElement(xfrm, _qn("a", "off"))
        off.set("x", "0")
        off.set("y", "0")
        ext = ET.SubElement(xfrm, _qn("a", "ext"))
        ext.set("cx", str(cx))
        ext.set("cy", str(cy))
        prst = ET.SubElement(sppr, _qn("a", "prstGeom"))
        prst.set("prst", "rect")
        ET.SubElement(prst, _qn("a", "avLst"))

        if caption:
            self.add_paragraph(
                caption,
                align="center",
                italic=True,
                size_half_points=20,  # 10pt
                spacing_before_twips=0,
                spacing_after_twips=240,
            )

    def _finalize_document_xml(self) -> bytes:
        # Append section properties as the last element in body.
        sect = ET.SubElement(self._body, _qn("w", "sectPr"))
        pg_sz = ET.SubElement(sect, _qn("w", "pgSz"))
        pg_sz.set(_qn("w", "w"), str(12240))  # 8.5 in * 1440 twips
        pg_sz.set(_qn("w", "h"), str(15840))  # 11 in * 1440 twips
        pg_mar = ET.SubElement(sect, _qn("w", "pgMar"))
        pg_mar.set(_qn("w", "top"), str(1440))
        pg_mar.set(_qn("w", "right"), str(1440))
        pg_mar.set(_qn("w", "bottom"), str(1440))
        pg_mar.set(_qn("w", "left"), str(1440))
        pg_mar.set(_qn("w", "header"), str(720))
        pg_mar.set(_qn("w", "footer"), str(720))
        pg_mar.set(_qn("w", "gutter"), str(0))

        doc = ET.Element(
            _qn("w", "document"),
            {
                "xmlns:w": NS["w"],
                "xmlns:r": NS["r"],
                "xmlns:wp": NS["wp"],
                "xmlns:a": NS["a"],
                "xmlns:pic": NS["pic"],
            },
        )
        doc.append(self._body)
        return ET.tostring(doc, encoding="utf-8", xml_declaration=True)

    def _document_rels_xml(self) -> bytes:
        rels = ET.Element(_qn("rels", "Relationships"), {"xmlns": NS["rels"]})
        for img in self._images:
            rel = ET.SubElement(rels, _qn("rels", "Relationship"))
            rel.set("Id", img.rid)
            rel.set("Type", "http://schemas.openxmlformats.org/officeDocument/2006/relationships/image")
            rel.set("Target", f"media/{img.filename}")
        return ET.tostring(rels, encoding="utf-8", xml_declaration=True)

    def write(self, output_path: Path) -> None:
        output_path.parent.mkdir(parents=True, exist_ok=True)

        # Root relationships.
        root_rels = ET.Element(_qn("rels", "Relationships"), {"xmlns": NS["rels"]})
        r1 = ET.SubElement(root_rels, _qn("rels", "Relationship"))
        r1.set("Id", "rId1")
        r1.set("Type", "http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument")
        r1.set("Target", "word/document.xml")
        r2 = ET.SubElement(root_rels, _qn("rels", "Relationship"))
        r2.set("Id", "rId2")
        r2.set("Type", "http://schemas.openxmlformats.org/package/2006/relationships/metadata/core-properties")
        r2.set("Target", "docProps/core.xml")
        r3 = ET.SubElement(root_rels, _qn("rels", "Relationship"))
        r3.set("Id", "rId3")
        r3.set("Type", "http://schemas.openxmlformats.org/officeDocument/2006/relationships/extended-properties")
        r3.set("Target", "docProps/app.xml")
        root_rels_xml = ET.tostring(root_rels, encoding="utf-8", xml_declaration=True)

        # Content types.
        types = ET.Element("Types", {"xmlns": "http://schemas.openxmlformats.org/package/2006/content-types"})
        ET.SubElement(
            types,
            "Default",
            {
                "Extension": "rels",
                "ContentType": "application/vnd.openxmlformats-package.relationships+xml",
            },
        )
        ET.SubElement(types, "Default", {"Extension": "xml", "ContentType": "application/xml"})
        ET.SubElement(types, "Default", {"Extension": "png", "ContentType": "image/png"})
        ET.SubElement(
            types,
            "Override",
            {
                "PartName": "/word/document.xml",
                "ContentType": "application/vnd.openxmlformats-officedocument.wordprocessingml.document.main+xml",
            },
        )
        ET.SubElement(
            types,
            "Override",
            {
                "PartName": "/docProps/core.xml",
                "ContentType": "application/vnd.openxmlformats-package.core-properties+xml",
            },
        )
        ET.SubElement(
            types,
            "Override",
            {
                "PartName": "/docProps/app.xml",
                "ContentType": "application/vnd.openxmlformats-officedocument.extended-properties+xml",
            },
        )
        content_types_xml = ET.tostring(types, encoding="utf-8", xml_declaration=True)

        # Core props.
        now = _dt.datetime.now(tz=_dt.timezone.utc).replace(microsecond=0)
        w3cdtf = now.isoformat().replace("+00:00", "Z")
        core = ET.Element(
            _qn("cp", "coreProperties"),
            {
                "xmlns:cp": NS["cp"],
                "xmlns:dc": NS["dc"],
                "xmlns:dcterms": NS["dcterms"],
                "xmlns:dcmitype": NS["dcmitype"],
                "xmlns:xsi": NS["xsi"],
            },
        )
        ET.SubElement(core, _qn("dc", "title")).text = self._title
        ET.SubElement(core, _qn("dc", "creator")).text = "Codex"
        ET.SubElement(core, _qn("cp", "lastModifiedBy")).text = "Codex"
        created = ET.SubElement(core, _qn("dcterms", "created"))
        created.set(_qn("xsi", "type"), "dcterms:W3CDTF")
        created.text = w3cdtf
        modified = ET.SubElement(core, _qn("dcterms", "modified"))
        modified.set(_qn("xsi", "type"), "dcterms:W3CDTF")
        modified.text = w3cdtf
        core_xml = ET.tostring(core, encoding="utf-8", xml_declaration=True)

        # App props.
        app = ET.Element(
            _qn("ep", "Properties"),
            {"xmlns": NS["ep"], "xmlns:vt": NS["vt"]},
        )
        ET.SubElement(app, _qn("ep", "Application")).text = "Codex"
        app_xml = ET.tostring(app, encoding="utf-8", xml_declaration=True)

        document_xml = self._finalize_document_xml()
        document_rels_xml = self._document_rels_xml()

        with zipfile.ZipFile(output_path, "w", compression=zipfile.ZIP_DEFLATED) as zf:
            zf.writestr("[Content_Types].xml", content_types_xml)
            zf.writestr("_rels/.rels", root_rels_xml)
            zf.writestr("docProps/core.xml", core_xml)
            zf.writestr("docProps/app.xml", app_xml)
            zf.writestr("word/document.xml", document_xml)
            zf.writestr("word/_rels/document.xml.rels", document_rels_xml)
            for img in self._images:
                zf.writestr(f"word/media/{img.filename}", img.data)


def _read_references(path: Path) -> list[str]:
    lines: list[str] = []
    for raw in path.read_text(encoding="utf-8").splitlines():
        s = raw.strip()
        if s:
            lines.append(s)
    return lines


def _iter_markdown_blocks(lines: Iterable[str]) -> Iterable[tuple[str, ...]]:
    # Blocks are: ("heading", level, text), ("paragraph", text), ("bullet", text),
    # ("image", caption, path), ("references",).
    para_parts: list[str] = []

    def flush_para():
        nonlocal para_parts
        if para_parts:
            text = " ".join(para_parts).strip()
            if text:
                yield ("paragraph", text)
            para_parts = []

    image_re = re.compile(r"^!\[(.*?)\]\((.*?)\)\s*$")

    for raw in lines:
        line = raw.rstrip("\n")
        stripped = line.strip()

        if stripped == "{{REFERENCES_TOP50}}":
            yield from flush_para()
            yield ("references",)
            continue

        mimg = image_re.match(stripped)
        if mimg:
            yield from flush_para()
            caption = mimg.group(1).strip()
            path = mimg.group(2).strip()
            yield ("image", caption, path)
            continue

        if stripped.startswith("#"):
            mh = re.match(r"^(#+)\s+(.*)$", stripped)
            if mh:
                yield from flush_para()
                level = len(mh.group(1))
                text = mh.group(2).strip()
                yield ("heading", str(level), text)
                continue

        if stripped.startswith("- "):
            yield from flush_para()
            yield ("bullet", stripped[2:].strip())
            continue

        if not stripped:
            yield from flush_para()
            continue

        para_parts.append(stripped)

    yield from flush_para()


def build_docx(draft_path: Path, references_path: Path, output_path: Path) -> None:
    refs = _read_references(references_path)
    blocks = list(_iter_markdown_blocks(draft_path.read_text(encoding="utf-8").splitlines()))
    title = None
    for b in blocks:
        if b and b[0] == "heading" and b[1] == "1":
            title = b[2]
            break
    if title is None:
        title = draft_path.stem

    w = _DocxWriter(title=title)

    for block in blocks:
        kind = block[0]
        if kind == "heading":
            level = int(block[1])
            text = block[2]
            if level == 1:
                w.add_paragraph(
                    text,
                    align="center",
                    bold=True,
                    size_half_points=56,  # 28pt
                    spacing_before_twips=0,
                    spacing_after_twips=360,
                )
            elif level == 2:
                w.add_paragraph(
                    text,
                    bold=True,
                    size_half_points=32,  # 16pt
                    spacing_before_twips=360,
                    spacing_after_twips=120,
                )
            else:
                w.add_paragraph(
                    text,
                    bold=True,
                    size_half_points=28,  # 14pt
                    spacing_before_twips=240,
                    spacing_after_twips=120,
                )
        elif kind == "paragraph":
            w.add_paragraph(block[1], size_half_points=22, spacing_after_twips=180)
        elif kind == "bullet":
            w.add_paragraph(
                f"- {block[1]}",
                size_half_points=22,
                indent_left_twips=360,
                spacing_after_twips=60,
            )
        elif kind == "image":
            caption, path_str = block[1], block[2]
            img_path = (draft_path.parent.parent / path_str).resolve() if not Path(path_str).is_absolute() else Path(path_str)
            w.add_image(img_path, caption=caption)
        elif kind == "references":
            for line in refs:
                w.add_paragraph(line, size_half_points=20, spacing_after_twips=60)
        else:
            raise ValueError(f"Unknown block: {block}")

    w.write(output_path)


def main(argv: list[str] | None = None) -> int:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description="Build manuscript .docx from manuscript/draft.md")
    parser.add_argument("--draft", type=Path, default=root / "manuscript" / "draft.md")
    parser.add_argument("--references", type=Path, default=root / "manuscript" / "references_top50.txt")
    parser.add_argument("--output", type=Path, default=root / "doravirine_resistance_mechanisms_draft.docx")
    args = parser.parse_args(argv)

    build_docx(args.draft, args.references, args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
