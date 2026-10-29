"""Edit a .docx in place without making Word ask to repair it.

Word is strict about how `word/document.xml` is serialised. Editing it with
`xml.etree.ElementTree` and writing it back naively changes three things that
Word notices:

1. `ET.tostring(..., xml_declaration=True)` emits
   ``<?xml version='1.0' encoding='UTF-8'?>`` -- single quotes and no
   ``standalone="yes"``. Word writes the double-quoted form with
   ``standalone`` and complains about the other.
2. Any namespace prefix not registered before serialisation is rewritten to
   ``ns0:``, ``ns1:`` and so on. A Word document declares about 35 of them.
3. The zip entries are rewritten with whatever compression and ordering the
   writer chooses; `[Content_Types].xml` in particular should come first.

`edit_document` handles all three. Use it instead of hand-rolling the
zip/ElementTree dance:

    from nnrti.docx_edit import edit_document

    with edit_document("paper.docx") as root:
        for p in root.iter(qn("w:p")):
            ...

The context manager writes the file only if the body exits without raising, so
a failed edit leaves the original untouched.
"""
from __future__ import annotations

import re
import shutil
import zipfile
from contextlib import contextmanager
from pathlib import Path
from xml.etree import ElementTree as ET

W = "http://schemas.openxmlformats.org/wordprocessingml/2006/main"

DECLARATION = b'<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\r\n'


def qn(tag: str) -> str:
    """'w:p' -> '{...wordprocessingml...}p', for ElementTree's tag syntax."""
    prefix, _, local = tag.partition(":")
    if not local:
        return tag
    uri = _NAMESPACES.get(prefix)
    if uri is None:
        raise KeyError(f"unknown namespace prefix {prefix!r}")
    return f"{{{uri}}}{local}"


_NAMESPACES: dict[str, str] = {"w": W}


def register_document_namespaces(xml: str) -> dict[str, str]:
    """Register every prefix the document declares, so none is rewritten to ns0."""
    # Four of them (a, a14, asvg, pic) are declared on inner drawing elements
    # rather than on <w:document>, so scanning the root tag alone leaves them
    # unregistered and ElementTree renames them to ns0, ns1, ...
    found = dict(re.findall(r'xmlns:([\w\d]+)="([^"]+)"', xml))
    for prefix, uri in found.items():
        ET.register_namespace(prefix, uri)
        _NAMESPACES[prefix] = uri
    return found


def serialise(root: ET.Element, original_root_tag: str | None = None) -> bytes:
    """Word-compatible bytes for document.xml.

    Two things have to survive serialisation.

    ElementTree drops namespace declarations it judges unused. That is valid
    XML and fatal here: ``mc:Ignorable`` on <w:document> names ten prefixes,
    and Word reports the file as unreadable if any of them is no longer
    declared. So every declaration from the original opening tag is merged
    back in, along with attributes such as ``mc:Ignorable`` itself.

    But the original tag cannot simply replace the generated one: four
    prefixes (a, a14, asvg, pic) are declared on inner drawing elements rather
    than on the root, and ElementTree hoists them to the root when it writes.
    Dropping its tag would leave those unbound. Merge, do not replace.
    """
    body = ET.tostring(root, encoding="utf-8", xml_declaration=False)
    if not original_root_tag:
        return DECLARATION + body

    generated = re.match(rb"<[^>]*>", body)
    if not generated:
        return DECLARATION + body
    gen_tag = generated.group(0).decode("utf-8")

    gen_ns = dict(re.findall(r'xmlns:([\w\d]+)="([^"]+)"', gen_tag))
    extra = []
    for prefix, uri in re.findall(r'xmlns:([\w\d]+)="([^"]+)"', original_root_tag):
        if prefix not in gen_ns:
            extra.append(f'xmlns:{prefix}="{uri}"')
    # non-namespace attributes of the root, mc:Ignorable above all
    for name, value in re.findall(r'(?<![\w:])((?!xmlns)[\w\d]+:[\w\d]+)="([^"]*)"',
                                  original_root_tag):
        if f"{name}=" not in gen_tag:
            extra.append(f'{name}="{value}"')

    if extra:
        gen_tag = gen_tag.rstrip(">").rstrip("/") + " " + " ".join(extra) + ">"
    return DECLARATION + gen_tag.encode("utf-8") + body[generated.end():]


@contextmanager
def edit_document(path: str | Path, part: str = "word/document.xml", backup: bool = True):
    """Yield the parsed root of `part`; write the archive back on clean exit."""
    path = Path(path)
    with zipfile.ZipFile(path) as zin:
        raw = zin.read(part).decode("utf-8")
        items = zin.infolist()
        payload = {item.filename: zin.read(item.filename) for item in items}

    register_document_namespaces(raw)
    match = re.search(r"<(?:\w+:)?document[^>]*>", raw)
    original_root_tag = match.group(0) if match else None
    root = ET.fromstring(raw)

    yield root                      # caller mutates the tree

    payload[part] = serialise(root, original_root_tag)
    if backup:
        shutil.copy2(path, path.with_suffix(path.suffix + ".bak"))

    tmp = path.with_suffix(path.suffix + ".tmp")
    with zipfile.ZipFile(tmp, "w", zipfile.ZIP_DEFLATED) as zout:
        # keep the original entry order; Word expects [Content_Types].xml first
        for item in items:
            zout.writestr(item, payload[item.filename])
    tmp.replace(path)


def check_document(path: str | Path) -> list[str]:
    """Return a list of problems Word would object to. Empty means clean."""
    path = Path(path)
    problems: list[str] = []
    with zipfile.ZipFile(path) as z:
        names = z.namelist()
        if "[Content_Types].xml" not in names:
            problems.append("missing [Content_Types].xml")
        elif names[0] != "[Content_Types].xml":
            problems.append(f"[Content_Types].xml is not first (found {names[0]})")
        raw = z.read("word/document.xml")

    if not raw.startswith(b'<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'):
        problems.append("XML declaration is not Word's form")
    text = raw.decode("utf-8", "replace")
    try:
        ET.fromstring(text)
    except ET.ParseError as exc:
        problems.append(f"document.xml does not parse: {exc}")
    if re.search(r"<ns\d+:", text):
        problems.append("generic ns0-style namespace prefixes present")
    match = re.search(r"<(?:\w+:)?document[^>]*>", text)
    if match:
        root_tag = match.group(0)
        if "mc:Ignorable" not in root_tag:
            problems.append("mc:Ignorable lost from the root element")
        else:
            declared = set(re.findall(r"xmlns:([\w\d]+)=", root_tag))
            ignorable = set(re.search(r'mc:Ignorable="([^"]*)"', root_tag).group(1).split())
            missing = sorted(ignorable - declared)
            if missing:
                problems.append(
                    "mc:Ignorable names undeclared prefixes (Word will offer to "
                    f"repair the file): {missing}"
                )

    starts = set(re.findall(r'<w:bookmarkStart[^>]*w:id="(\d+)"', text))
    ends = set(re.findall(r'<w:bookmarkEnd[^>]*w:id="(\d+)"', text))
    if starts ^ ends:
        problems.append(f"unmatched bookmarks: {sorted(starts ^ ends)[:6]}")

    begins = text.count('w:fldCharType="begin"')
    closes = text.count('w:fldCharType="end"')
    if begins != closes:
        problems.append(f"unbalanced fields: {begins} begin, {closes} end")

    return problems
