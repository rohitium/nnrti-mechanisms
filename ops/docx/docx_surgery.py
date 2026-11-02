#!/usr/bin/env python3
"""Run-aware surgery on a .docx paragraph, safe around Zotero fields.

Word splits a sentence across arbitrarily many runs, often mid-word, so text
has to be located in the paragraph's reconstructed visible text and then mapped
back onto the runs the match actually touches. Citation fields make that
delicate: a Zotero citation is five runs (begin, instrText, separate, the
visible superscript, end) and only the fourth carries text, so a span that
includes a citation number necessarily spans the whole field.

This module therefore does not merely tolerate fields, it moves them. Deleting
a span returns the field element groups it contained, so a sentence that is
being relocated to another paragraph can take its citations with it, intact and
still bound to the same Zotero item. That is the operation the reorganisation
needs and the one a naive text rewrite gets wrong.

Every helper works on lxml/ElementTree elements of word/document.xml. Nothing
here writes a file; the caller owns packing.
"""
from __future__ import annotations

import copy
import re
import shutil
import zipfile
from pathlib import Path
import xml.etree.ElementTree as ET

W = "{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"
XML_SPACE = "{http://www.w3.org/XML/1998/namespace}space"
ET.register_namespace("w", W[1:-1])
DOC = "word/document.xml"
RED = "FF0000"


# --------------------------------------------------------------------------- io
def load(path: Path) -> tuple[ET.ElementTree, bytes]:
    with zipfile.ZipFile(path) as z:
        return ET.ElementTree(ET.fromstring(z.read(DOC))), z.read(DOC)


def save(src: Path, dst: Path, tree: ET.ElementTree) -> None:
    """Repack, replacing only word/document.xml and preserving everything else."""
    blob = ET.tostring(tree.getroot(), encoding="UTF-8", xml_declaration=True)
    if src != dst:
        shutil.copy2(src, dst)
    tmp = dst.with_suffix(dst.suffix + ".tmp")
    with zipfile.ZipFile(dst) as zin, zipfile.ZipFile(tmp, "w", zipfile.ZIP_DEFLATED) as zout:
        for item in zin.infolist():
            zout.writestr(item, blob if item.filename == DOC else zin.read(item.filename))
    tmp.replace(dst)


def body(tree: ET.ElementTree) -> ET.Element:
    return tree.getroot().find(W + "body")


def paragraphs(tree: ET.ElementTree) -> list[ET.Element]:
    return list(body(tree).iter(W + "p"))


# ------------------------------------------------------------------- inspection
def run_text(run: ET.Element) -> str:
    return "".join(t.text or "" for t in run.iter(W + "t"))


def para_text_raw(par: ET.Element) -> str:
    """Visible text exactly as the runs hold it, double spaces and all.

    Spans passed to the editing helpers must be sliced from this, not from
    para_text: Word leaves stray double spaces behind, and a span normalised
    for reading will not be found in the runs.
    """
    return "".join(run_text(r) for r in par.iter(W + "r"))


def para_text(par: ET.Element) -> str:
    """Visible text, citation superscripts included, as a reader sees it."""
    return re.sub(r"\s+", " ", para_text_raw(par)).strip()


def _runs(par: ET.Element) -> list[ET.Element]:
    """Direct-child runs. Deeper nesting would break offset mapping, so refuse it."""
    direct = [c for c in par if c.tag == W + "r"]
    if len(direct) != len(list(par.iter(W + "r"))):
        raise RuntimeError("paragraph has runs nested below the top level")
    return direct


def field_count(par: ET.Element) -> tuple[int, int]:
    xml = ET.tostring(par, encoding="unicode")
    return xml.count('w:fldCharType="begin"'), xml.count('w:fldCharType="end"')


def _locate(par: ET.Element, needle: str) -> tuple[list[ET.Element], int, int, int, int]:
    runs = _runs(par)
    texts = [run_text(r) for r in runs]
    joined = "".join(texts)
    hits = [m.start() for m in re.finditer(re.escape(needle), joined)]
    if len(hits) != 1:
        raise LookupError(f"{len(hits)} matches for {needle[:60]!r}")
    start, end = hits[0], hits[0] + len(needle)
    pos = 0
    si = so = ei = eo = None
    for i, t in enumerate(texts):
        nxt = pos + len(t)
        if si is None and pos <= start < nxt:
            si, so = i, start - pos
        if pos < end <= nxt:
            ei, eo = i, end - pos
        pos = nxt
    if si is None or ei is None:
        raise LookupError("could not map the match onto runs")
    return runs, si, so, ei, eo


def _set_text(run: ET.Element, text: str) -> None:
    ts = list(run.iter(W + "t"))
    if not ts:
        if text:
            raise RuntimeError("cannot put text in a run that has no w:t")
        return
    ts[0].text = text
    ts[0].set(XML_SPACE, "preserve")
    for extra in ts[1:]:
        extra.text = ""


def make_run(text: str, template: ET.Element | None = None, red: bool = True) -> ET.Element:
    """A run carrying `text`, red by default so every addition is visible as one."""
    run = ET.Element(W + "r")
    if red:
        rpr = ET.SubElement(run, W + "rPr")
        ET.SubElement(rpr, W + "color").set(W + "val", RED)
    elif template is not None:
        src = template.find(W + "rPr")
        if src is not None:
            run.append(copy.deepcopy(src))
    t = ET.SubElement(run, W + "t")
    t.text = text
    t.set(XML_SPACE, "preserve")
    return run


def redden_field(group: list[ET.Element]) -> list[ET.Element]:
    """Colour a field's visible superscript red, leaving its instrText alone."""
    for run in group:
        if run.find(W + "t") is not None:
            rpr = run.find(W + "rPr")
            if rpr is None:
                rpr = ET.Element(W + "rPr")
                run.insert(0, rpr)
            for old in rpr.findall(W + "color"):
                rpr.remove(old)
            ET.SubElement(rpr, W + "color").set(W + "val", RED)
    return group


def new_citation(group: list[ET.Element], citation_id: str, shown: str) -> list[ET.Element]:
    """Clone a field as a fresh citation of the same item(s).

    Zotero keys each field by citationID, so a clone needs a new one or the two
    instances are the same citation as far as a refresh is concerned.
    """
    out = [copy.deepcopy(r) for r in group]
    for run in out:
        for instr in run.iter(W + "instrText"):
            if instr.text and "citationID" in instr.text:
                instr.text = re.sub(r'"citationID":"[^"]*"',
                                    f'"citationID":"{citation_id}"', instr.text)
        if run.find(W + "t") is not None:
            _set_text(run, shown)
    return redden_field(out)


# ---------------------------------------------------------------- modifications
def delete_span(par: ET.Element, needle: str) -> list[list[ET.Element]]:
    """Remove `needle`; return the field groups it contained, ready to re-use."""
    runs, si, so, ei, eo = _locate(par, needle)
    doomed = runs[si + 1:ei] if ei > si else []
    fields, group, depth = [], [], 0
    for run in doomed:
        xml = ET.tostring(run, encoding="unicode")
        if 'w:fldCharType="begin"' in xml:
            depth += 1
        if depth:
            group.append(run)
        if 'w:fldCharType="end"' in xml:
            depth -= 1
            if depth == 0:
                fields.append(group)
                group = []
    if depth:
        raise RuntimeError("the span ends inside a field; refusing")

    if si == ei:
        t = run_text(runs[si])
        _set_text(runs[si], t[:so] + t[eo:])
    else:
        _set_text(runs[si], run_text(runs[si])[:so])
        _set_text(runs[ei], run_text(runs[ei])[eo:])
        for run in doomed:
            par.remove(run)
    return fields


def replace_span(par: ET.Element, old: str, new: str, red: bool = True) -> None:
    """Swap `old` for `new`, which arrives as its own red run."""
    runs, si, so, ei, eo = _locate(par, old)
    inside = runs[si + 1:ei] if ei > si else []
    if any('w:fldChar' in ET.tostring(r, encoding="unicode") for r in inside):
        raise RuntimeError("replacement span crosses a field; delete and rebuild instead")
    if si == ei:
        t = run_text(runs[si])
        head, tail = t[:so], t[eo:]
        _set_text(runs[si], head)
        at = list(par).index(runs[si])
        par.insert(at + 1, make_run(new, runs[si], red))
        if tail:
            par.insert(at + 2, make_run(tail, runs[si], red=False))
    else:
        _set_text(runs[si], run_text(runs[si])[:so])
        _set_text(runs[ei], run_text(runs[ei])[eo:])
        for run in inside:
            par.remove(run)
        par.insert(list(par).index(runs[si]) + 1, make_run(new, runs[si], red))


def append_run(par: ET.Element, text: str, red: bool = True) -> None:
    par.append(make_run(text, None, red))


def build_paragraph(parts, style_from: ET.Element) -> ET.Element:
    """A new paragraph in the same style, from ('text', s) and ('field', group)."""
    par = ET.Element(W + "p")
    ppr = style_from.find(W + "pPr")
    if ppr is not None:
        par.append(copy.deepcopy(ppr))
    for kind, value in parts:
        if kind == "text":
            par.append(make_run(value, None, red=True))
        elif kind == "field":
            for run in value:
                par.append(run)
        else:
            raise ValueError(kind)
    return par


def insert_after(container: ET.Element, anchor: ET.Element, new_paras) -> None:
    at = list(container).index(anchor)
    for offset, par in enumerate(new_paras, start=1):
        container.insert(at + offset, par)


def find_para(tree: ET.ElementTree, anchor: str) -> ET.Element:
    """The one paragraph whose visible text contains `anchor`.

    Paragraph indices are not a usable address: whether table cells count
    changes the numbering, and an insertion renumbers everything after it. An
    anchor that matches exactly one paragraph is stable across both.
    """
    hits = [p for p in paragraphs(tree) if anchor in para_text(p)]
    if len(hits) != 1:
        raise LookupError(f"{len(hits)} paragraphs contain {anchor[:70]!r}")
    return hits[0]


def parent_of(tree: ET.ElementTree, par: ET.Element) -> ET.Element:
    for candidate in tree.getroot().iter():
        if par in list(candidate):
            return candidate
    raise LookupError("paragraph has no parent")


def find_paras(tree: ET.ElementTree, anchor: str) -> list[ET.Element]:
    """Every paragraph containing `anchor`, in document order.

    The Response quotes one Discussion paragraph twice, under two different
    reviewer comments, so some anchors legitimately match more than once and
    the two copies have to be addressed by position.
    """
    return [p for p in paragraphs(tree) if anchor in para_text(p)]
