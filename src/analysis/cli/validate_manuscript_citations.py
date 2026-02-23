#!/usr/bin/env python3
"""Validate numeric citations and report investigator coverage.

Checks:
  - Parses in-text citations like [12] and ranges like [12-15].
  - Confirms references section contains 1..N with no gaps.
  - Confirms every reference number is cited at least once.
  - Reports counts of references containing: Arnold, Hughes, Stammers.
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path


def _expand_citation_token(token: str) -> list[int]:
    t = token.strip()
    if not t:
        return []
    if "-" in t:
        a, b = t.split("-", 1)
        try:
            start = int(a.strip())
            end = int(b.strip())
        except ValueError:
            return []
        if start <= 0 or end <= 0:
            return []
        if end < start:
            start, end = end, start
        return list(range(start, end + 1))
    try:
        n = int(t)
    except ValueError:
        return []
    return [n] if n > 0 else []


def parse_intext_citations(text: str) -> set[int]:
    # Capture bracket groups, then split on commas/semicolons.
    # Examples:
    #   [1], [1-3], [1,4,7-9]
    cite_re = re.compile(r"\[([0-9,\s;\-]+)\]")
    out: set[int] = set()
    for m in cite_re.finditer(text):
        inner = m.group(1)
        for part in re.split(r"[,\s;]+", inner.strip()):
            out.update(_expand_citation_token(part))
    return {n for n in out if n > 0}


def parse_reference_numbers(lines: list[str]) -> dict[int, str]:
    # Matches "12. ..." at line start.
    ref_re = re.compile(r"^\s*(\d+)\.\s+(.*\S)\s*$")
    refs: dict[int, str] = {}
    for line in lines:
        m = ref_re.match(line)
        if not m:
            continue
        n = int(m.group(1))
        refs[n] = m.group(2).strip()
    return refs


def main() -> int:
    parser = argparse.ArgumentParser(description="Validate manuscript numeric citations + investigator coverage.")
    parser.add_argument(
        "--manuscript",
        type=Path,
        default=Path("manuscript/doravirine_resistance_mechanisms_draft.md"),
    )
    args = parser.parse_args()

    if not args.manuscript.exists():
        print(f"ERROR: missing manuscript: {args.manuscript}", file=sys.stderr)
        return 2

    text = args.manuscript.read_text(encoding="utf-8")
    lines = text.splitlines()

    cited = parse_intext_citations(text)

    refs = parse_reference_numbers(lines)
    if not refs:
        print("ERROR: no references found (expected lines like '1. ...')", file=sys.stderr)
        return 2

    max_ref = max(refs)
    expected = set(range(1, max_ref + 1))
    missing_numbers = sorted(expected - set(refs))
    if missing_numbers:
        print(f"ERROR: reference list has gaps; missing numbers: {missing_numbers}", file=sys.stderr)
        return 2

    uncited = sorted(expected - cited)
    if uncited:
        print(f"ERROR: {len(uncited)} references are never cited in text. First 25: {uncited[:25]}", file=sys.stderr)
        return 2

    # Investigator coverage.
    def _count_name(name: str) -> int:
        patt = re.compile(rf"\b{name}\b", flags=re.IGNORECASE)
        return sum(1 for _n, ref in refs.items() if patt.search(ref))

    arnold = _count_name("Arnold")
    hughes = _count_name("Hughes")
    stammers = _count_name("Stammers")

    print(f"OK: citations={len(cited)} unique cited numbers; references={max_ref} entries; all cited.")
    print(f"Coverage: Arnold={arnold} refs; Hughes={hughes} refs; Stammers={stammers} refs.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

