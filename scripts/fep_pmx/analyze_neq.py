#!/usr/bin/env python3
"""Run pmx BAR analysis on completed NEQ switch work values."""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.fep_pmx.config import FEP_PMX_ROOT, NEQ_TEMPERATURE_K


def _collect_dgdl(neq_dir: Path, *, direction: str) -> list[Path]:
    pattern = "fwd" if direction == "fwd" else "rev"
    files = sorted((neq_dir / "switches").glob(f"{pattern}_*/dgdl.xvg"))
    if not files:
        raise FileNotFoundError(f"No {pattern} dgdl.xvg files under {neq_dir / 'switches'}")
    return files


def analyze_neq_leg(
    leg_id: str,
    *,
    phase: str,
    replicate: int = 1,
    temperature_k: float = NEQ_TEMPERATURE_K,
    output_dir: Path | None = None,
) -> Path:
    neq_dir = FEP_PMX_ROOT / "legs" / leg_id / phase / f"rep_{replicate:02d}" / "neq"
    fwd = _collect_dgdl(neq_dir, direction="fwd")
    rev = _collect_dgdl(neq_dir, direction="rev")

    out_dir = output_dir or (neq_dir / "analysis")
    out_dir.mkdir(parents=True, exist_ok=True)
    summary = out_dir / "bar_summary.txt"

    cmd = [
        "pmx",
        "analyse",
        "-fA",
        *[str(p) for p in fwd],
        "-fB",
        *[str(p) for p in rev],
        "-t",
        str(temperature_k),
        "-m",
        "bar",
        "-o",
        str(out_dir / "bar"),
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
    summary.write_text(
        "\n".join(
            [
                f"leg_id={leg_id}",
                f"phase={phase}",
                f"replicate={replicate}",
                f"n_fwd={len(fwd)}",
                f"n_rev={len(rev)}",
                f"returncode={proc.returncode}",
                "",
                "stdout:",
                proc.stdout,
                "",
                "stderr:",
                proc.stderr,
            ]
        )
        + "\n"
    )
    if proc.returncode != 0:
        raise RuntimeError(f"pmx analyse failed; see {summary}")

    meta = {
        "leg_id": leg_id,
        "phase": phase,
        "replicate": replicate,
        "temperature_k": temperature_k,
        "n_fwd": len(fwd),
        "n_rev": len(rev),
        "summary": str(summary),
        "bar_prefix": str(out_dir / "bar"),
    }
    (out_dir / "analysis.json").write_text(json.dumps(meta, indent=2) + "\n")
    return summary


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Analyze pmx NEQ BAR results for one leg.")
    parser.add_argument("--leg", required=True)
    parser.add_argument("--phase", choices=("holo", "apo"), required=True)
    parser.add_argument("--replicate", type=int, default=1)
    parser.add_argument("--temperature-k", type=float, default=NEQ_TEMPERATURE_K)
    args = parser.parse_args(argv)

    try:
        summary = analyze_neq_leg(
            args.leg,
            phase=args.phase,
            replicate=args.replicate,
            temperature_k=args.temperature_k,
        )
    except (FileNotFoundError, RuntimeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    print(f"Wrote analysis summary: {summary}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
