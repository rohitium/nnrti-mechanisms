#!/usr/bin/env python3
from pathlib import Path
import argparse
import csv


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Rewrite absolute local paths in fep_manifest.csv to a Sherlock root"
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        default=Path("results/fep_manifest.csv"),
        help="Manifest CSV to rewrite (default: results/fep_manifest.csv)",
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
        help="New Sherlock repo root prefix",
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
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
