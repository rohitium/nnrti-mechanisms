#!/usr/bin/env python3
"""Effect sizes and confidence intervals for every structural feature.

Reviewer 2 (comment 3) asks for effect sizes with confidence intervals on the
comparisons behind the mechanistic claims. Rather than treat the two
comparisons named in the comment as special cases, this applies one procedure
to every genotype and every structural feature, in both the bound and
unliganded states, so a reader can check any claim in the paper the same way.

The unit of replication is the replicate, not the frame. Frames within a
replicate are autocorrelated and treating them as independent would shrink
every interval by an order of magnitude; three replicates means three numbers.

For each genotype and feature, against wild type:

    shift      difference of replicate means, in the feature's own units
    sem        standard error of that difference
    t, df, p   Welch's unequal-variance t test (df ~ 4 for 3 vs 3)
    ci95       shift +/- t(0.975, df) x sem

A note on thresholds, because it is easy to get wrong at n = 3. A shift of
twice its standard error is NOT significant here: with df = 4 that is p = 0.12.
The 95% interval excludes zero only beyond 2.78 standard errors. Both are
reported so a claim can be checked against whichever the text uses.

Usage
-----
    PYTHONPATH=src python -m nnrti.cli.compute_structural_effect_sizes
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

#: Columns that identify a row rather than measure something.
_ID_COLUMNS = {"structure_type", "genotype", "replicate", "dor_fold_reduction"}

#: Dihedrals, which must be averaged and differenced on the circle.
_CIRCULAR_SUFFIX = ("chi1_Tyr181_deg", "chi1_Tyr188_deg")


def _circ_mean(deg: np.ndarray) -> float:
    r = np.deg2rad(np.asarray(deg, float))
    r = r[np.isfinite(r)]
    if r.size == 0:
        return np.nan
    return float(np.rad2deg(np.arctan2(np.sin(r).mean(), np.cos(r).mean())))


def _centre_circular(values: np.ndarray, centre: float) -> np.ndarray:
    """Re-express angles as signed deviations from a centre, so that the
    ordinary t machinery applies to a quantity that is no longer wrapped."""
    return ((np.asarray(values, float) - centre + 180.0) % 360.0) - 180.0


def compare(mut: np.ndarray, wt: np.ndarray, circular: bool) -> dict:
    """Welch comparison of two small samples of replicate means."""
    from scipy import stats

    mut = np.asarray(mut, float)
    wt = np.asarray(wt, float)
    mut, wt = mut[np.isfinite(mut)], wt[np.isfinite(wt)]
    if mut.size < 2 or wt.size < 2:
        return {}

    if circular:
        # Centre both samples on the wild-type circular mean before testing:
        # a dihedral near +/-180 would otherwise have a meaningless variance.
        c = _circ_mean(wt)
        mut, wt = _centre_circular(mut, c), _centre_circular(wt, c)

    t, p = stats.ttest_ind(mut, wt, equal_var=False)
    vm, vw = mut.var(ddof=1), wt.var(ddof=1)
    nm, nw = mut.size, wt.size
    sem = float(np.sqrt(vm / nm + vw / nw))
    # Welch-Satterthwaite degrees of freedom.
    df = float((vm / nm + vw / nw) ** 2 /
               ((vm / nm) ** 2 / (nm - 1) + (vw / nw) ** 2 / (nw - 1))) \
        if sem > 0 else np.nan
    shift = float(mut.mean() - wt.mean())
    crit = float(stats.t.ppf(0.975, df)) if np.isfinite(df) else np.nan

    return {
        "n_mut": nm, "n_wt": nw, "shift": shift, "sem": sem,
        "t": float(t), "df": df, "p": float(p),
        "ci95_low": shift - crit * sem, "ci95_high": shift + crit * sem,
        "sem_ratio": abs(shift) / sem if sem > 0 else np.nan,
    }


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--per-replicate", type=Path,
                   default=Path("results/analysis/structural_metrics/tables/"
                                "structural_metrics_per_replicate.csv"))
    p.add_argument("--output-dir", type=Path,
                   default=Path("results/analysis/structural_metrics/tables"))
    args = p.parse_args(argv)

    d = pd.read_csv(args.per_replicate)
    features = [c for c in d.columns if c not in _ID_COLUMNS]

    rows = []
    for state, g in d.groupby("structure_type"):
        wt = g[g.genotype == "WT"]
        if wt.empty:
            continue
        for feat in features:
            if g[feat].isna().all():
                continue
            circular = feat in _CIRCULAR_SUFFIX
            for gt, gg in g[g.genotype != "WT"].groupby("genotype"):
                stats_ = compare(gg[feat].to_numpy(), wt[feat].to_numpy(), circular)
                if not stats_:
                    continue
                rows.append({"structure_type": state, "genotype": gt,
                             "feature": feat, **stats_})

    out = pd.DataFrame(rows)
    for c in ("shift", "sem", "t", "df", "ci95_low", "ci95_high", "sem_ratio"):
        out[c] = out[c].round(4)
    out["p"] = out["p"].round(5)
    out["significant_95"] = (out.p < 0.05)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    f = args.output_dir / "structural_effect_sizes.csv"
    out.to_csv(f, index=False)

    n = len(out)
    print(f"{n} genotype x feature comparisons, both states")
    print(f"  |shift| > 2.00 x SEM (p ~ 0.12 at df=4): "
          f"{int((out.sem_ratio > 2.0).sum()):4d}  ({100*(out.sem_ratio>2).mean():.0f}%)")
    print(f"  95% CI excludes zero (p < 0.05)        : "
          f"{int(out.significant_95.sum()):4d}  ({100*out.significant_95.mean():.0f}%)")
    print()
    for state, g in out.groupby("structure_type"):
        print(f"  {state}: {int(g.significant_95.sum())}/{len(g)} significant at 95%")
    print(f"\nwrote {f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
