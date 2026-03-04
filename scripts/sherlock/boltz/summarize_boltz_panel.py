#!/usr/bin/env python3
"""Summarize Boltz affinity replicate JSON files for a mutation panel."""

import argparse
import csv
import glob
import json
import math
import statistics as st
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, List, Tuple

TCRIT_95 = {
    1: 12.706,
    2: 4.303,
    3: 3.182,
    4: 2.776,
    5: 2.571,
    6: 2.447,
    7: 2.365,
    8: 2.306,
    9: 2.262,
    10: 2.228,
    11: 2.201,
    12: 2.179,
    13: 2.160,
    14: 2.145,
    15: 2.131,
    16: 2.120,
    17: 2.110,
    18: 2.101,
    19: 2.093,
    20: 2.086,
    21: 2.080,
    22: 2.074,
    23: 2.069,
    24: 2.064,
    25: 2.060,
    26: 2.056,
    27: 2.052,
    28: 2.048,
    29: 2.045,
    30: 2.042,
}


def tcrit95(n):
    if n <= 1:
        return float("nan")
    return TCRIT_95.get(n - 1, 1.96)


def summarize(values):
    n = len(values)
    mean = st.mean(values)
    if n <= 1:
        return mean, 0.0, 0.0
    sd = st.stdev(values)
    ci = tcrit95(n) * sd / math.sqrt(n)
    return mean, sd, ci


def mutation_from_path(path):
    # Expected: .../<mutation>/replicates/affinity_seedXXXX.json
    try:
        return path.parent.parent.name
    except Exception:
        return "unknown"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--glob",
        dest="input_glob",
        default="results/boltz/control_panel/*/replicates/affinity_seed*.json",
        help="Glob pattern for replicate affinity JSON files",
    )
    parser.add_argument(
        "--out-csv",
        type=Path,
        default=None,
        help="Optional output CSV path for summary rows",
    )
    args = parser.parse_args()

    files = sorted(glob.glob(args.input_glob))
    if not files:
        raise FileNotFoundError(f"No files found for pattern: {args.input_glob}")

    by_mut = defaultdict(
        lambda: {"log10_uM": [], "prob": [], "uM": [], "nM": []}
    )  # type: Dict[str, Dict[str, List[float]]]

    for file_path in files:
        path = Path(file_path)
        mutation = mutation_from_path(path)
        with path.open() as fh:
            payload = json.load(fh)

        log10_um = float(payload["affinity_pred_value"])
        prob = float(payload["affinity_probability_binary"])
        um = 10 ** log10_um
        nm = um * 1000.0

        by_mut[mutation]["log10_uM"].append(log10_um)
        by_mut[mutation]["prob"].append(prob)
        by_mut[mutation]["uM"].append(um)
        by_mut[mutation]["nM"].append(nm)

    rows = []  # type: List[Dict[str, Any]]
    for mutation in sorted(by_mut):
        x = by_mut[mutation]
        n = len(x["log10_uM"])
        m_log, sd_log, ci_log = summarize(x["log10_uM"])
        m_prob, sd_prob, ci_prob = summarize(x["prob"])
        m_um, sd_um, ci_um = summarize(x["uM"])
        m_nm, sd_nm, ci_nm = summarize(x["nM"])
        rows.append(
            {
                "mutation": mutation,
                "n": n,
                "mean_log10_uM": m_log,
                "sd_log10_uM": sd_log,
                "ci95_log10_uM": ci_log,
                "mean_prob": m_prob,
                "sd_prob": sd_prob,
                "ci95_prob": ci_prob,
                "mean_ic50_uM": m_um,
                "sd_ic50_uM": sd_um,
                "ci95_ic50_uM": ci_um,
                "mean_ic50_nM": m_nm,
                "sd_ic50_nM": sd_nm,
                "ci95_ic50_nM": ci_nm,
            }
        )

    headers = [
        "mutation",
        "n",
        "mean_log10_uM",
        "sd_log10_uM",
        "ci95_log10_uM",
        "mean_prob",
        "sd_prob",
        "ci95_prob",
        "mean_ic50_uM",
        "sd_ic50_uM",
        "ci95_ic50_uM",
        "mean_ic50_nM",
        "sd_ic50_nM",
        "ci95_ic50_nM",
    ]

    print("\t".join(headers))
    for row in rows:
        print(
            "{}\t{}\t{:.6g}\t{:.6g}\t{:.6g}\t{:.6g}\t{:.6g}\t{:.6g}\t{:.6g}\t{:.6g}\t{:.6g}\t{:.6g}\t{:.6g}\t{:.6g}".format(
                row["mutation"],
                row["n"],
                row["mean_log10_uM"],
                row["sd_log10_uM"],
                row["ci95_log10_uM"],
                row["mean_prob"],
                row["sd_prob"],
                row["ci95_prob"],
                row["mean_ic50_uM"],
                row["sd_ic50_uM"],
                row["ci95_ic50_uM"],
                row["mean_ic50_nM"],
                row["sd_ic50_nM"],
                row["ci95_ic50_nM"],
            )
        )

    if args.out_csv is not None:
        args.out_csv.parent.mkdir(parents=True, exist_ok=True)
        with args.out_csv.open("w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=headers)
            writer.writeheader()
            writer.writerows(rows)
        print(f"\nWrote summary CSV: {args.out_csv}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
