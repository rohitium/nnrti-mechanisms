#!/usr/bin/env python3
from __future__ import annotations

import os
import sys
from pathlib import Path


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def main(argv: list[str] | None = None) -> int:
    import argparse

    root = _repo_root()
    os.environ.setdefault("MPLCONFIGDIR", str(root / ".mplconfig"))
    sys.path.insert(0, str(root))

    import matplotlib.pyplot as plt

    from src.susceptibility_io import load_dor_susceptibilities

    parser = argparse.ArgumentParser(description="Plot ranked doravirine fold-reduction from the susceptibility workbook")
    parser.add_argument(
        "--input",
        type=Path,
        default=root / "data" / "DRM-susceptibilities.csv.xlsx",
        help="Path to susceptibility workbook (.xlsx).",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=root / "results" / "plots" / "dor_fold_reduction_ranked.png",
        help="Output PNG path.",
    )
    args = parser.parse_args(argv)

    df = load_dor_susceptibilities(args.input)
    df = df.sort_values("dor_fold_reduction", ascending=False).reset_index(drop=True)

    mutations = df["mutation"].tolist()
    values = df["dor_fold_reduction"].astype(float).tolist()

    # Dynamic height so labels stay readable for small/large panels.
    fig_h = max(3.5, 0.45 * len(mutations) + 1.2)
    fig_w = 9.0
    plt.figure(figsize=(fig_w, fig_h))

    ax = plt.gca()
    y = list(range(len(mutations)))
    ax.barh(y, values, color="#2b6cb0", alpha=0.9)
    ax.set_yticks(y, labels=mutations)
    ax.invert_yaxis()

    ax.set_xlabel("DOR median fold reduction")
    ax.set_title("Doravirine (DOR) susceptibility landscape for selected DRMs")

    xmax = max(values) if values else 1.0
    xpad = 0.03 * xmax
    for yi, v in zip(y, values):
        label = f"{v:.1f}".rstrip("0").rstrip(".")
        ax.text(v + xpad, yi, label, va="center", ha="left", fontsize=10)

    ax.set_xlim(0.0, xmax * 1.15)
    ax.grid(axis="x", linestyle=":", linewidth=0.8, alpha=0.5)
    plt.tight_layout()

    args.output.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(args.output, dpi=300, bbox_inches="tight")
    plt.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

