#!/usr/bin/env python3
"""Run pmx BAR/CGI/Jarzynski on completed NEQ switch work values.

Wraps ``pmx analyse`` for one leg/phase/replicate: dumps per-switch integrated
work values (``integ_fwd.dat`` / ``integ_rev.dat``), pmx's own work-distribution
plot, the full ``results.txt``, and a machine-readable ``analysis.json`` with the
parsed BAR/CGI/Jarzynski free energies (in kcal/mol).

``combine_neq.py`` and ``qc_neq.py`` build on the ``analysis.json`` this writes.
"""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.fep_pmx.config import FEP_PMX_ROOT, NEQ_TEMPERATURE_K

OUTPUT_UNITS = "kcal/mol"
KJ_PER_KCAL = 4.184


def _neq_dir(leg_id: str, phase: str, replicate: int) -> Path:
    return FEP_PMX_ROOT / "legs" / leg_id / phase / f"rep_{replicate:02d}" / "neq"


def _collect_dgdl(neq_dir: Path, *, direction: str) -> list[Path]:
    files = sorted((neq_dir / "switches").glob(f"{direction}_*/dgdl.xvg"))
    if not files:
        raise FileNotFoundError(f"No {direction} dgdl.xvg files under {neq_dir / 'switches'}")
    return [p.resolve() for p in files]


def _find_float(pattern: str, text: str) -> float | None:
    """First capture group of ``pattern`` as float; None if absent or non-finite."""
    match = re.search(pattern, text, flags=re.MULTILINE)
    if not match:
        return None
    try:
        value = float(match.group(1))
    except ValueError:
        return None
    if value != value or value in (float("inf"), float("-inf")):  # nan / inf
        return None
    return value


def parse_pmx_results(results_text: str) -> dict[str, float | None]:
    """Parse the estimator free energies from a ``pmx analyse`` results file.

    Matches the exact ``_tee`` format strings in pmx ``analyze_dhdl.py`` (values
    already in the requested output units, here kcal/mol).
    """
    num = r"(-?[\d.]+)"
    return {
        "bar_dg": _find_float(rf"^\s*BAR:\s*dG\s*=\s*{num}", results_text),
        "bar_err_analytical": _find_float(
            rf"^\s*BAR:\s*Std Err \(analytical\)\s*=\s*{num}", results_text
        ),
        "bar_err_boot": _find_float(rf"^\s*BAR:\s*Std Err \(bootstrap\)\s*=\s*{num}", results_text),
        "bar_conv": _find_float(rf"^\s*BAR:\s*Conv\s*=\s*{num}", results_text),
        "cgi_dg": _find_float(rf"^\s*CGI:\s*dG\s*=\s*{num}", results_text),
        "jarz_dg_mean": _find_float(rf"^\s*JARZ:\s*dG Mean\s*=\s*{num}", results_text),
    }


def analyze_neq_leg(
    leg_id: str,
    *,
    phase: str,
    replicate: int = 1,
    temperature_k: float = NEQ_TEMPERATURE_K,
    nboots: int = 100,
    output_dir: Path | None = None,
    force: bool = False,
) -> dict:
    """Run pmx analyse for one leg/phase/replicate and return the parsed result."""
    neq_dir = _neq_dir(leg_id, phase, replicate)
    out_dir = output_dir or (neq_dir / "analysis")
    analysis_json = out_dir / "analysis.json"
    if analysis_json.is_file() and not force:
        return json.loads(analysis_json.read_text())

    fwd = _collect_dgdl(neq_dir, direction="fwd")
    rev = _collect_dgdl(neq_dir, direction="rev")

    out_dir.mkdir(parents=True, exist_ok=True)
    results_txt = out_dir / "results.txt"
    integ_fwd = out_dir / "integ_fwd.dat"
    integ_rev = out_dir / "integ_rev.dat"
    work_plot = out_dir / "work_dist.png"

    # Absolute paths throughout so the invocation is cwd-independent.
    cmd = [
        "pmx", "analyse",
        "-fA", *[str(p) for p in fwd],
        "-fB", *[str(p) for p in rev],
        "-m", "cgi", "bar", "jarz",
        "-t", str(temperature_k),
        "-b", str(nboots),
        "--units", "kcal",
        "-o", str(results_txt.resolve()),
        "-oA", str(integ_fwd.resolve()),
        "-oB", str(integ_rev.resolve()),
        "-w", str(work_plot.resolve()),
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
    (out_dir / "pmx_stdout.txt").write_text(proc.stdout + "\n---stderr---\n" + proc.stderr + "\n")
    if proc.returncode != 0:
        raise RuntimeError(f"pmx analyse failed (rc={proc.returncode}); see {out_dir / 'pmx_stdout.txt'}")

    parsed = parse_pmx_results(results_txt.read_text() if results_txt.is_file() else proc.stdout)
    meta = {
        "leg_id": leg_id,
        "phase": phase,
        "replicate": replicate,
        "temperature_k": temperature_k,
        "units": OUTPUT_UNITS,
        "n_fwd": len(fwd),
        "n_rev": len(rev),
        "nboots": nboots,
        "results_txt": str(results_txt),
        "integ_fwd": str(integ_fwd),
        "integ_rev": str(integ_rev),
        "work_plot": str(work_plot),
        **parsed,
    }
    analysis_json.write_text(json.dumps(meta, indent=2) + "\n")
    return meta


def ensure_leg_analysis(
    leg_id: str,
    *,
    phase: str,
    replicate: int,
    temperature_k: float = NEQ_TEMPERATURE_K,
    nboots: int = 100,
    auto: bool = True,
) -> dict:
    """Return the parsed analysis dict, running pmx analyse first if needed."""
    analysis_json = _neq_dir(leg_id, phase, replicate) / "analysis" / "analysis.json"
    if analysis_json.is_file():
        return json.loads(analysis_json.read_text())
    if not auto:
        raise FileNotFoundError(
            f"Missing {analysis_json}; run analyze_neq.py for {leg_id} {phase} rep{replicate}"
        )
    return analyze_neq_leg(
        leg_id, phase=phase, replicate=replicate,
        temperature_k=temperature_k, nboots=nboots,
    )


def read_work_values_kcal(integ_dat: Path) -> list[float]:
    """Read per-switch integrated work values (pmx dumps them in kJ/mol)."""
    works: list[float] = []
    for line in Path(integ_dat).read_text().splitlines():
        line = line.strip()
        if not line or line.startswith(("#", "@")):
            continue
        try:
            works.append(float(line.split()[-1]) / KJ_PER_KCAL)
        except (ValueError, IndexError):
            continue
    return works


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Analyze pmx NEQ results for one leg.")
    parser.add_argument("--leg", required=True)
    parser.add_argument("--phase", choices=("holo", "apo"), required=True)
    parser.add_argument("--replicate", type=int, default=1)
    parser.add_argument("--temperature-k", type=float, default=NEQ_TEMPERATURE_K)
    parser.add_argument("--nboots", type=int, default=100)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args(argv)

    try:
        meta = analyze_neq_leg(
            args.leg, phase=args.phase, replicate=args.replicate,
            temperature_k=args.temperature_k, nboots=args.nboots, force=args.force,
        )
    except (FileNotFoundError, RuntimeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    print(
        f"{args.leg} {args.phase} rep{args.replicate}: "
        f"BAR dG = {meta['bar_dg']} ± {meta['bar_err_analytical']} {meta['units']} "
        f"(n_fwd={meta['n_fwd']}, n_rev={meta['n_rev']})"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
