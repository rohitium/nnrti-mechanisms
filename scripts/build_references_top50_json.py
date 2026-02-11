#!/usr/bin/env python3
from __future__ import annotations

import json
import re
from pathlib import Path


def _tags_for(citation: str) -> list[str]:
    c = citation.lower()
    tags: list[str] = []
    if "doravirine" in c or "mk-1439" in c:
        tags.append("dor")
    if "resistance" in c or "susceptib" in c or "mutat" in c:
        tags.append("resistance")
    if "reverse transcriptase" in c or "rt-" in c:
        tags.append("rt")
    if "structure" in c or "crystal" in c or "nmr" in c:
        tags.append("structure")
    if "molecular dynamics" in c or "openmm" in c or "mdanalysis" in c or "mm/pbsa" in c or "mm/gbsa" in c:
        tags.append("methods")
    # De-duplicate while preserving order.
    out: list[str] = []
    seen = set()
    for t in tags:
        if t not in seen:
            seen.add(t)
            out.append(t)
    return out


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    txt_path = root / "manuscript" / "references_top50.txt"
    json_path = root / "manuscript" / "references_top50.json"

    if not txt_path.exists():
        raise FileNotFoundError(f"Missing: {txt_path}")

    records: list[dict] = []
    for raw in txt_path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line:
            continue
        m = re.match(r"^(\d+)\.\s+(.*)$", line)
        if not m:
            raise ValueError(f"Unrecognized reference line: {line!r}")
        ref_id = int(m.group(1))
        citation = m.group(2).strip()

        doi = None
        mdoi = re.search(r"doi:([^\s]+)", citation, flags=re.IGNORECASE)
        if mdoi:
            doi = mdoi.group(1).rstrip(".").strip()

        pmid = None
        mpmid = re.search(r"PMID:(\d+)", citation, flags=re.IGNORECASE)
        if mpmid:
            pmid = mpmid.group(1)

        year = None
        myear = re.search(r"\.\s*(\d{4})[;\.]", citation)
        if myear:
            year = int(myear.group(1))

        records.append(
            {
                "id": ref_id,
                "citation": citation,
                "year": year,
                "doi": doi,
                "pmid": pmid,
                "tags": _tags_for(citation),
                # Citation counts were not collected deterministically in this offline build.
                "citation_count": None,
            }
        )

    if len(records) != 50:
        raise ValueError(f"Expected 50 references, found {len(records)}")
    if sorted(r["id"] for r in records) != list(range(1, 51)):
        raise ValueError("Reference ids must be exactly 1..50")

    json_path.write_text(json.dumps(records, indent=2, sort_keys=False), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
