#!/usr/bin/env python3
"""Free-energy profile G(λ) along the alchemical morphing coordinate.

For each leg/phase, reconstruct free energy as a function of λ (0 = WT, 1 = mutant)
from the per-switch dH/dλ trajectories (switches/*/dgdl.xvg) — the alchemical
analog of the PMF-along-path in Serra et al. 2025 (JCTC 21, 2079), their Fig. 2/3.

Per switch, the cumulative work is  W(λ) = ∫₀^λ (dH/dλ') dλ'  (λ = t / t_switch).
From the forward switches we plot, per leg/phase:
  * ⟨W_f(λ)⟩  — mean forward cumulative work (upper envelope)
  * G_f(λ)    — forward Jarzynski free energy, −kT ln⟨exp(−W_f(λ)/kT)⟩
  * dissipation(λ) = ⟨W_f(λ)⟩ − G_f(λ)  (shaded) — where irreversibility builds up
The converged BAR ΔG (from analysis.json) is drawn at λ = 1 as the reference.

Run where the dgdl.xvg live (Sherlock). Outputs per leg:
  results/analysis/fep_pmx/lambda_profiles/<leg>.png  and  <leg>.csv

  python3 scripts/fep_pmx/plot_lambda_profile.py --legs wt_to_V106A wt_to_K103N
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.fep_pmx.config import FEP_PMX_ROOT, NEQ_TEMPERATURE_K

KB_KCAL = 0.0019872041  # kcal/(mol*K)
KJ_PER_KCAL = 4.184


def _work_profile(dgdl: Path, grid: np.ndarray) -> np.ndarray:
    """Cumulative work W(λ) [kcal/mol] on `grid` for one switch (forward frame).

    Reads (time, dH/dλ); λ = normalized time in [0, 1]; cumulative trapezoid of
    dH/dλ over λ, converted kJ->kcal, interpolated onto the common grid.
    """
    d = np.loadtxt(dgdl, comments=["#", "@"])
    t = d[:, 0]
    dhdl = d[:, -1]
    lam = (t - t[0]) / (t[-1] - t[0])
    w = np.concatenate([[0.0], np.cumsum(0.5 * (dhdl[1:] + dhdl[:-1]) * np.diff(lam))])
    return np.interp(grid, lam, w) / KJ_PER_KCAL


def _jarzynski(work_by_switch: np.ndarray, kt: float) -> np.ndarray:
    """G(λ) = -kT ln mean_i exp(-W_i(λ)/kT), log-sum-exp for stability. Axis 0 = switch."""
    a = -work_by_switch / kt
    m = a.max(axis=0)
    return -kt * (m + np.log(np.mean(np.exp(a - m), axis=0)))


def profile_leg(leg_id: str, *, replicates: range, npts: int, out_dir: Path,
                temperature_k: float) -> Path | None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    kt = KB_KCAL * temperature_k
    grid = np.linspace(0.0, 1.0, npts)
    phases = ("holo", "apo")
    data: dict[str, dict] = {}
    for phase in phases:
        wf_all = []
        bar_dgs = []
        for rep in replicates:
            neq = FEP_PMX_ROOT / "legs" / leg_id / phase / f"rep_{rep:02d}" / "neq"
            fwd = sorted((neq / "switches").glob("fwd_*/dgdl.xvg"))
            for f in fwd:
                try:
                    wf_all.append(_work_profile(f, grid))
                except Exception as exc:  # noqa: BLE001
                    print(f"  warn {f}: {exc}")
            aj = neq / "analysis" / "analysis.json"
            if aj.is_file():
                bd = json.loads(aj.read_text()).get("bar_dg")
                if bd is not None:
                    bar_dgs.append(bd)
        if not wf_all:
            continue
        wf = np.vstack(wf_all)
        wf_mean = wf.mean(axis=0)
        g_fwd = _jarzynski(wf, kt)
        data[phase] = {
            "wf_mean": wf_mean, "g_fwd": g_fwd, "diss": wf_mean - g_fwd,
            "n_switch": wf.shape[0], "bar_dg": float(np.mean(bar_dgs)) if bar_dgs else None,
        }
        print(f"{leg_id} {phase}: n_switch={wf.shape[0]}  "
              f"⟨W_f(1)⟩={wf_mean[-1]:+.1f}  G_f(1)[Jarz]={g_fwd[-1]:+.1f}  "
              f"BAR ΔG={data[phase]['bar_dg']}  dissipation(1)={data[phase]['diss'][-1]:.1f} kcal/mol")
    if not data:
        return None

    ncols = len(data)
    fig, axes = plt.subplots(1, ncols, figsize=(6.4 * ncols, 4.6), squeeze=False)
    for c, (phase, d) in enumerate(data.items()):
        ax = axes[0][c]
        ax.plot(grid, d["wf_mean"], color="#2c6fbb", lw=1.8, label=r"$\langle W_f(\lambda)\rangle$ (forward work)")
        ax.plot(grid, d["g_fwd"], color="#1a7a3a", lw=1.8, label=r"$G_f(\lambda)$ (forward Jarzynski)")
        ax.fill_between(grid, d["g_fwd"], d["wf_mean"], color="#d1642f", alpha=0.18,
                        label="dissipation")
        if d["bar_dg"] is not None:
            ax.plot(1.0, d["bar_dg"], "k*", ms=12, zorder=5, label=f"BAR ΔG = {d['bar_dg']:.2f}")
        ax.axhline(0, color="0.7", lw=0.8)
        ax.set_xlabel(r"$\lambda$  (0 = WT  →  1 = mutant)")
        ax.set_ylabel("free energy / work (kcal/mol)")
        ax.set_title(f"{leg_id}  ·  {phase}  ({d['n_switch']} fwd switches)", fontsize=10)
        ax.legend(fontsize=8, loc="best")
    fig.suptitle(f"{leg_id}: free-energy profile along λ", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    out_dir.mkdir(parents=True, exist_ok=True)
    png = out_dir / f"{leg_id}.png"
    fig.savefig(png, dpi=200)
    plt.close(fig)

    # light CSV for re-plotting / verification
    csv = out_dir / f"{leg_id}.csv"
    cols = ["lambda"]
    arrs = [grid]
    for phase, d in data.items():
        cols += [f"{phase}_Wf_mean", f"{phase}_G_fwd", f"{phase}_dissipation"]
        arrs += [d["wf_mean"], d["g_fwd"], d["diss"]]
    np.savetxt(csv, np.column_stack(arrs), header=",".join(cols), delimiter=",", comments="")
    return png


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description="G(λ) free-energy profile per leg.")
    p.add_argument("--legs", nargs="+", required=True)
    p.add_argument("--replicates", type=int, default=3)
    p.add_argument("--npoints", type=int, default=101)
    p.add_argument("--temperature-k", type=float, default=NEQ_TEMPERATURE_K)
    p.add_argument("--output-dir", type=Path, default=FEP_PMX_ROOT / "lambda_profiles")
    args = p.parse_args(argv)
    n = 0
    for leg in args.legs:
        out = profile_leg(leg, replicates=range(1, args.replicates + 1), npts=args.npoints,
                          out_dir=args.output_dir, temperature_k=args.temperature_k)
        if out:
            print(f"wrote {out}\n")
            n += 1
        else:
            print(f"skip {leg}: no fwd dgdl found\n")
    print(f"{n} leg profile(s) written to {args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
