#!/usr/bin/env python3
"""Effect sizes with confidence intervals for the key structural comparisons.

Reviewer 2 (comment 3) noted that several mechanistic claims rest on
differences reported only as means with replicate SEMs, and asked for effect
sizes with confidence intervals for specific comparisons -- V106A against
V106I/M, and G190A against G190E.

The unit of replication is the trajectory, not the frame. Frames within a
replicate are strongly autocorrelated, so treating them as independent would
shrink the interval by roughly the square root of the frame count and state a
precision the data does not have. Every statistic here is therefore computed
from the three replicate means per genotype.

With n = 3 per group the interval, not the point estimate, is the honest
output. Two are reported:

  * Welch's t interval on the difference of replicate means. Primary, because
    it does not assume equal variance and is defined for n = 3.
  * A percentile bootstrap resampling replicates within each genotype.
    Secondary and deliberately labelled crude: three values admit only 27
    distinct resamples, so the interval is granular.

Hedges' g is reported with the small-sample correction, which matters at this n.

Usage
-----
    PYTHONPATH=src python -m nnrti.cli.compute_effect_sizes
"""

from __future__ import annotations

import argparse
import itertools
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

#: comparisons the reviewer named, plus WT references for context
COMPARISONS = [
    ("V106A", "V106I"),
    ("V106A", "V106M"),
    ("G190A", "G190E"),
    ("wt", "V106A"),
    ("wt", "G190E"),
    ("wt", "Y188L"),
]

#: per-frame metrics from mechanism_coordinates.csv
FRAME_METRICS = {
    "dor_to_S105_mindist": "DOR-Ser105 minimum distance (A)",
    "chl_ring_burial": "Chlorocyanophenyl ring burial (contacts)",
    "res103_bb_to_triazN": "Lys103 C=O to triazolinone N (A)",
    "dor_to_res106_mindist": "DOR-residue 106 minimum distance (A)",
}


def _hedges_g(a: np.ndarray, b: np.ndarray) -> float:
    na, nb = len(a), len(b)
    sp = np.sqrt(((na - 1) * a.var(ddof=1) + (nb - 1) * b.var(ddof=1)) / (na + nb - 2))
    if sp == 0:
        return float("nan")
    d = (b.mean() - a.mean()) / sp
    # Hedges' small-sample correction; at n=3 per group this is ~0.85
    j = 1.0 - 3.0 / (4.0 * (na + nb) - 9.0)
    return float(d * j)


def _bootstrap_ci(a: np.ndarray, b: np.ndarray, n_boot: int, rng) -> tuple[float, float]:
    diffs = np.empty(n_boot)
    for i in range(n_boot):
        diffs[i] = rng.choice(b, len(b), replace=True).mean() - rng.choice(a, len(a), replace=True).mean()
    return float(np.percentile(diffs, 2.5)), float(np.percentile(diffs, 97.5))


def compare(a: np.ndarray, b: np.ndarray, n_boot: int, rng) -> dict:
    """Difference b - a, with Welch and bootstrap intervals."""
    diff = float(b.mean() - a.mean())
    t = stats.ttest_ind(b, a, equal_var=False)
    # Welch interval on the difference
    se = np.sqrt(b.var(ddof=1) / len(b) + a.var(ddof=1) / len(a))
    df = (b.var(ddof=1) / len(b) + a.var(ddof=1) / len(a)) ** 2 / (
        (b.var(ddof=1) / len(b)) ** 2 / (len(b) - 1) + (a.var(ddof=1) / len(a)) ** 2 / (len(a) - 1)
    ) if se > 0 else float("nan")
    tcrit = stats.t.ppf(0.975, df) if np.isfinite(df) and df > 0 else float("nan")
    lo, hi = (diff - tcrit * se, diff + tcrit * se) if np.isfinite(tcrit) else (float("nan"),) * 2
    blo, bhi = _bootstrap_ci(a, b, n_boot, rng)
    return {
        "mean_a": float(a.mean()), "mean_b": float(b.mean()),
        "n_a": len(a), "n_b": len(b),
        "difference": diff,
        "welch_ci_low": lo, "welch_ci_high": hi,
        "welch_p": float(t.pvalue), "welch_df": float(df) if np.isfinite(df) else float("nan"),
        "boot_ci_low": blo, "boot_ci_high": bhi,
        "hedges_g": _hedges_g(a, b),
        "crosses_zero": bool(not np.isfinite(lo) or (lo <= 0 <= hi)),
    }


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--mechanisms", type=Path,
                   default=Path("results/analysis/mechanisms/mechanism_coordinates.csv"))
    p.add_argument("--moieties", type=Path,
                   default=Path("results/analysis/mechanisms/dor_moiety_contacts_per_replicate.csv"))
    p.add_argument("--pocket", type=Path,
                   default=Path("results/analysis/pocket_volume/tables/pocket_volume_per_rep.csv"))
    p.add_argument("--output-dir", type=Path, default=Path("results/analysis/effect_sizes"))
    p.add_argument("--n-boot", type=int, default=20000)
    p.add_argument("--seed", type=int, default=20261020)
    args = p.parse_args(argv)

    rng = np.random.default_rng(args.seed)
    out = args.output_dir / "tables"
    out.mkdir(parents=True, exist_ok=True)

    # Replicate-level means for each source.
    frames = pd.read_csv(args.mechanisms)
    rep_frame = frames.groupby(["safe_label", "replicate"])[list(FRAME_METRICS)].mean().reset_index()

    moi = pd.read_csv(args.moieties)
    moi["safe_label"] = moi["mutation"].str.replace("+", "_", regex=False).str.replace(" ", "", regex=False)
    moi.loc[moi.mutation.str.upper() == "WT", "safe_label"] = "wt"
    rep_moi = moi[["safe_label", "replicate", "pyridinone", "chlorocyanophenyl"]]

    pock = pd.read_csv(args.pocket)
    pock["safe_label"] = pock["mutation"].str.replace("+", "_", regex=False).str.replace(" ", "", regex=False)
    pock.loc[pock.mutation.str.upper() == "WT", "safe_label"] = "wt"
    rep_pock = pock[["safe_label", "replicate", "pocket_volume_mean"]]

    table = (rep_frame.merge(rep_moi, on=["safe_label", "replicate"], how="outer")
                      .merge(rep_pock, on=["safe_label", "replicate"], how="outer"))

    metrics = {**FRAME_METRICS,
               "pyridinone": "Pyridinone moiety contacts",
               "chlorocyanophenyl": "Chlorocyanophenyl moiety contacts",
               "pocket_volume_mean": "NNIBP volume (A^3)"}

    rows = []
    for ga, gb in COMPARISONS:
        for col, label in metrics.items():
            a = table.loc[table.safe_label == ga, col].dropna().to_numpy(float)
            b = table.loc[table.safe_label == gb, col].dropna().to_numpy(float)
            if len(a) < 2 or len(b) < 2:
                continue
            r = {"comparison": f"{ga} -> {gb}", "reference": ga, "genotype": gb,
                 "metric": col, "metric_label": label}
            r.update(compare(a, b, args.n_boot, rng))
            rows.append(r)

    res = pd.DataFrame(rows)
    path = out / "effect_sizes.csv"
    res.to_csv(path, index=False)

    print(f"{len(res)} comparisons -> {path}\n")
    show = res[["comparison", "metric", "difference", "welch_ci_low", "welch_ci_high",
                "hedges_g", "crosses_zero"]].copy()
    for c in ("difference", "welch_ci_low", "welch_ci_high", "hedges_g"):
        show[c] = show[c].round(2)
    print(show.to_string(index=False))
    n_sig = int((~res.crosses_zero).sum())
    print(f"\n{n_sig} of {len(res)} comparisons have a 95% Welch interval excluding zero.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
