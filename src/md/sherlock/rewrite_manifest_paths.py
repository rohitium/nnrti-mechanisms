#!/usr/bin/env python3
from pathlib import Path
import argparse
import csv
import json


def _rewrite_json_results(results_dir: Path, from_root: str, to_root: str) -> int:
    """Rewrite absolute paths inside result JSON files."""
    count = 0
    for jf in results_dir.rglob("*.json"):
        text = jf.read_text()
        if from_root in text:
            jf.write_text(text.replace(from_root, to_root))
            count += 1
    return count


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Rewrite absolute paths in manifest CSV (and optionally result JSONs)"
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        default=Path("results/md_manifest.csv"),
        help="Manifest CSV to rewrite (default: results/md_manifest.csv)",
    )
    parser.add_argument(
        "--from-root",
        type=str,
        default="/Users/rohitpro/Career/00_Github/nnrti-mechanisms",
        help="Old absolute root prefix to replace",
    )
    parser.add_argument(
        "--to-root",
        type=Path,
        default=Path("/scratch/users/rsatija/nnrti-mechanisms"),
        help="New repo root prefix",
    )
    parser.add_argument(
        "--rewrite-jsons",
        type=Path,
        default=None,
        help="Also rewrite JSON result files under this directory",
    )
    args = parser.parse_args()

    manifest = args.manifest
    if not manifest.exists():
        raise FileNotFoundError(f"Manifest not found: {manifest}")

    rows = []
    with manifest.open("r", newline="") as handle:
        reader = csv.DictReader(handle)
        fieldnames = reader.fieldnames
        for row in reader:
            for key, value in row.items():
                if isinstance(value, str) and value.startswith(args.from_root):
                    row[key] = value.replace(args.from_root, str(args.to_root), 1)
            rows.append(row)

    with manifest.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f"Rewrote {len(rows)} rows in {manifest}")

    if args.rewrite_jsons and args.rewrite_jsons.is_dir():
        n = _rewrite_json_results(args.rewrite_jsons, args.from_root, str(args.to_root))
        if n:
            print(f"Rewrote {n} JSON result files in {args.rewrite_jsons}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
