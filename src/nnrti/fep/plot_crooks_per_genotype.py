#!/usr/bin/env python3
"""Per-genotype Crooks work-overlap figures (one PNG per mutation).

Splits the monolithic panel_crooks_overlap.png into a compact figure per
genotype: its leg(s) x phase(s), replicates pooled, forward vs reverse work with
the BAR dG and the overlap coefficient. Reuses qc_neq.qc_unit so the overlap
definition matches panel_qc.csv exactly. Uses cached analysis (no dgdl needed);
pass --force to re-run pmx analyse on Sherlock.

  python -m nnrti.fep.plot_crooks_per_genotype --targets K103N V106M ...
Outputs: results/analysis/fep_pmx/crooks_overlap/<genotype>.png
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[3]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from nnrti.fep.mutations import MANUSCRIPT_PLANS
from nnrti.fep.combine_neq import DEFAULT_TARGETS, EXPERIMENTAL_CSV, load_experimental
from nnrti.fep.config import FEP_PMX_ROOT, NEQ_TEMPERATURE_K
from nnrti.fep.qc_neq import qc_unit


def _sanitize(name: str) -> str:
    return name.replace("+", "_").replace("/", "_").replace(" ", "")


def plot_genotype(genotype: str, *, replicates: range, temperature_k: float, nboots: int,
                  auto: bool, force: bool, fold: float | None, out_dir: Path) -> Path | None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    legs = MANUSCRIPT_PLANS[genotype].legs
    phases = ("holo", "apo")
    # collect pooled work per (leg, phase)
    cells: dict[tuple[str, str], dict] = {}
    for leg in legs:
        for phase in phases:
            wf, wr, ov, dgs = [], [], [], []
            for rep in replicates:
                try:
                    u = qc_unit(leg.leg_id, phase, rep, temperature_k=temperature_k,
                                nboots=nboots, auto=auto, force=force)
                except Exception:
                    continue
                if u.get("_wf") is None or len(u["_wf"]) == 0:
                    continue
                wf.append(u["_wf"]); wr.append(u["_wr"]); ov.append(u["overlap"])
                if u.get("bar_dg") is not None:
                    dgs.append(u["bar_dg"])
            if wf:
                cells[(leg.leg_id, phase)] = {
                    "wf": np.concatenate(wf), "wr": np.concatenate(wr),
                    "overlap": float(np.mean(ov)), "bar_dg": float(np.mean(dgs)) if dgs else None,
                    "n": len(wf),
                }
    if not cells:
        return None

    nrows, ncols = len(legs), len(phases)
    fig, axes = plt.subplots(nrows, ncols, figsize=(5.4 * ncols, 3.2 * nrows), squeeze=False)
    for r, leg in enumerate(legs):
        for c, phase in enumerate(phases):
            ax = axes[r][c]
            cell = cells.get((leg.leg_id, phase))
            if cell is None:
                ax.set_visible(False)
                continue
            ax.hist(cell["wf"], bins=25, alpha=0.55, density=True, color="#2c6fbb", label=r"forward $W_f$")
            ax.hist(cell["wr"], bins=25, alpha=0.55, density=True, color="#d1642f", label=r"reverse $W_r$")
            if cell["bar_dg"] is not None:
                ax.axvline(cell["bar_dg"], color="0.2", lw=1.3, ls="--",
                           label=f"BAR ΔG = {cell['bar_dg']:.2f}")
            ax.set_title(f"{leg.leg_id}  ·  {phase}  (overlap {cell['overlap']:.2f}, {cell['n']} rep)",
                         fontsize=9)
            ax.set_xlabel("work (kcal/mol)"); ax.set_ylabel("density")
            ax.legend(fontsize=7)
    fold_s = f"fold {fold:g}" if fold is not None else "no experimental fold"
    fig.suptitle(f"{genotype}  —  Crooks work overlap  ({fold_s})", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / f"{_sanitize(genotype)}.png"
    fig.savefig(out, dpi=200)
    plt.close(fig)
    return out


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description="Per-genotype Crooks overlap figures.")
    p.add_argument("--targets", nargs="+", default=list(DEFAULT_TARGETS))
    p.add_argument("--replicates", type=int, default=3)
    p.add_argument("--temperature-k", type=float, default=NEQ_TEMPERATURE_K)
    p.add_argument("--nboots", type=int, default=100)
    p.add_argument("--no-auto", action="store_true")
    p.add_argument("--force", action="store_true")
    p.add_argument("--experimental-csv", type=Path, default=EXPERIMENTAL_CSV)
    p.add_argument("--output-dir", type=Path, default=FEP_PMX_ROOT / "crooks_overlap")
    args = p.parse_args(argv)

    fold = load_experimental(args.experimental_csv)
    reps = range(1, args.replicates + 1)
    written = 0
    for g in args.targets:
        if g not in MANUSCRIPT_PLANS:
            print(f"skip {g}: not in MANUSCRIPT_PLANS")
            continue
        out = plot_genotype(g, replicates=reps, temperature_k=args.temperature_k,
                            nboots=args.nboots, auto=not args.no_auto, force=args.force,
                            fold=fold.get(g), out_dir=args.output_dir)
        if out:
            print(f"wrote {out}")
            written += 1
        else:
            print(f"skip {g}: no work data found")
    print(f"\n{written} genotype figure(s) written to {args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
