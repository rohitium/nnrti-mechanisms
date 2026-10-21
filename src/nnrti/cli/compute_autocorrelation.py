#!/usr/bin/env python3
"""Autocorrelation times and effective sample sizes for the structural observables.

Reviewer 2 (comment 1) asked for convergence diagnostics beyond the block
averaging already computed, specifically autocorrelation analysis, and for a
statement of how limited sampling may bias the results.

Frames within a replicate are not independent. The integrated autocorrelation
time tau_int measures how many frames one effectively waits between
independent samples, and the effective sample size is

    ESS = N / (2 * tau_int)

tau_int is estimated with Geyer's initial positive sequence: sum the
autocorrelation over successive lag pairs and stop when a pair sums to a
negative value, which is the standard way to avoid the variance blow-up of
summing a noisy tail.

A caveat that belongs with the numbers: tau_int characterises motion *within*
the sampled basin. A replicate that never leaves its starting basin can have a
short autocorrelation time and a large ESS while still not sampling the
conformational landscape. These numbers bound the precision of the reported
means; they do not establish that the ensemble is the right one.

MM/GBSA and FEP are handled by their protocols rather than here. MM/GBSA
snapshots are drawn evenly across the post-equilibration region rather than
consecutively, and the FEP work values come from switches seeded from
separately drawn frames, so neither is a serially correlated time series.

Usage
-----
    PYTHONPATH=src python -m nnrti.cli.compute_autocorrelation
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

METRICS = [
    "dor_to_S105_mindist",
    "chl_ring_burial",
    "res103_bb_to_triazN",
    "dor_to_res106_mindist",
    "res179_to_dor_mindist",
    "y188_interplanar_deg",
]


def tau_int(x: np.ndarray) -> float:
    """Integrated autocorrelation time, Geyer initial positive sequence."""
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    n = len(x)
    if n < 8 or np.allclose(x, x[0]):
        return float("nan")
    x = x - x.mean()
    var = np.dot(x, x) / n
    if var == 0:
        return float("nan")
    # autocorrelation via FFT
    f = np.fft.rfft(x, 2 * n)
    acf = np.fft.irfft(f * np.conjugate(f))[:n] / (n * var)
    # Geyer: pair successive lags, stop when the pair turns negative
    total = 0.0
    for k in range(1, n - 1, 2):
        pair = acf[k] + acf[k + 1]
        if pair <= 0:
            break
        total += pair
    return float(1.0 + 2.0 * total)


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--mechanisms", type=Path,
                   default=Path("results/analysis/mechanisms/mechanism_coordinates.csv"))
    p.add_argument("--output-dir", type=Path, default=Path("results/analysis/md_convergence"))
    args = p.parse_args(argv)

    d = pd.read_csv(args.mechanisms)
    out = args.output_dir / "tables"
    out.mkdir(parents=True, exist_ok=True)

    rows = []
    for (lbl, rep), g in d.groupby(["safe_label", "replicate"]):
        g = g.sort_values("frame")
        dt = float(np.median(np.diff(g.time_ns))) if len(g) > 1 else float("nan")
        for m in METRICS:
            if m not in g:
                continue
            t = tau_int(g[m].to_numpy())
            n = int(g[m].notna().sum())
            rows.append({
                "safe_label": lbl, "replicate": rep, "metric": m,
                "n_frames": n, "frame_dt_ns": dt,
                "tau_int_frames": t,
                "tau_int_ns": t * dt if np.isfinite(t) else float("nan"),
                "ess": n / (2 * t) if np.isfinite(t) and t > 0 else float("nan"),
            })

    res = pd.DataFrame(rows)
    path = out / "autocorrelation.csv"
    res.to_csv(path, index=False)

    summ = (res.groupby("metric")
               .agg(tau_frames_median=("tau_int_frames", "median"),
                    tau_ns_median=("tau_int_ns", "median"),
                    ess_median=("ess", "median"),
                    ess_min=("ess", "min"),
                    n_series=("ess", "size"))
               .round(2))
    summ.to_csv(out / "autocorrelation_summary.csv")

    print(f"{len(res)} series -> {path}\n")
    print(summ.to_string())
    print("\nPer-genotype median ESS (across metrics and replicates):")
    print(res.groupby("safe_label").ess.median().round(1).sort_values().to_string())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
