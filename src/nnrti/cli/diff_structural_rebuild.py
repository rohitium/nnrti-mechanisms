#!/usr/bin/env python3
"""What a rebuild of the structural tables changes, cell by cell.

CLAUDE.md asks that a re-run reproduce the published numbers before its new
results are reported. For most of this panel it does, exactly -- so the useful
output is not "did anything change" but "which cells changed, by how much, and
does any sentence in the manuscript quote one of them".

Compares a rebuilt structural_metrics_per_replicate.csv against the backed-up
published copy, at three levels:

  replicate   every genotype x state x replicate x metric
  genotype    the three-replicate means and SEMs the tables print
  figure 5    the induced-fit double difference and its interval, which is
              what changes a figure cell from resolved to unresolved

Dihedrals are differenced on the circle at every level; an ordinary difference
of two angles near the wrap reports 300 degrees for a 60 degree move.

Usage
-----
    PYTHONPATH=src python -m nnrti.cli.diff_structural_rebuild \
        --published archive/analysis-updates/table-snapshot-YYYYMMDD/\
structural_metrics/tables/structural_metrics_per_replicate.csv
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from nnrti.cli.plot_structural_heatmap import CIRCULAR, DOR_METRICS, RT_METRICS

#: t(0.975) for the Welch df of each comparison, as plot_structural_heatmap uses.
T_INDUCED_FIT = 2.306


def circular_mean(deg) -> float:
    r = np.deg2rad(np.asarray(deg, dtype=float))
    return float(np.rad2deg(np.arctan2(np.sin(r).mean(), np.cos(r).mean())))


def circular_delta(a: float, b: float) -> float:
    return float((a - b + 180.0) % 360.0 - 180.0)


def circular_sem(vals) -> float:
    v = np.asarray(vals, dtype=float)
    v = v[np.isfinite(v)]
    if v.size < 2:
        return float("nan")
    m = circular_mean(v)
    dev = np.array([circular_delta(x, m) for x in v])
    return float(np.std(dev, ddof=1) / np.sqrt(dev.size))


def _agg(df: pd.DataFrame, col: str, circular: bool):
    g = df.groupby(["state", "genotype"])[col]
    if circular:
        return g.apply(circular_mean), g.apply(circular_sem)
    return g.mean(), g.apply(lambda s: float(s.std(ddof=1) / np.sqrt(s.notna().sum()))
                             if s.notna().sum() > 1 else np.nan)


def induced_fit(df: pd.DataFrame, col: str, circular: bool) -> pd.DataFrame:
    """(holo - apo)_mut - (holo - apo)_WT, with the error Figure 5 quotes."""
    mean, sem = _agg(df, col, circular)
    rows = []
    for geno in sorted({g for _s, g in mean.index}):
        if geno == "WT":
            continue
        keys = [("holo", geno), ("apo", geno), ("holo", "WT"), ("apo", "WT")]
        if any(k not in mean.index or not np.isfinite(mean[k]) for k in keys):
            continue
        if circular:
            d_mut = circular_delta(mean[keys[0]], mean[keys[1]])
            d_wt = circular_delta(mean[keys[2]], mean[keys[3]])
            value = circular_delta(d_mut, d_wt)
        else:
            value = (mean[keys[0]] - mean[keys[1]]) - (mean[keys[2]] - mean[keys[3]])
        err = float(np.sqrt(np.nansum([sem.get(k, np.nan) ** 2 for k in keys])))
        rows.append(dict(genotype=geno, value=value, err=err,
                         resolved=bool(abs(value) > T_INDUCED_FIT * err)))
    if not rows:
        # A DOR-contact metric exists only in the bound state, so it has no
        # induced-fit double difference at all. That is not a failure.
        return pd.DataFrame(columns=["genotype", "value", "err", "resolved"]
                            ).set_index("genotype")
    return pd.DataFrame(rows).set_index("genotype")


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--published", type=Path, required=True)
    p.add_argument("--rebuilt", type=Path,
                   default=Path("results/analysis/structural_metrics/tables/"
                                "structural_metrics_per_replicate.csv"))
    p.add_argument("--tol", type=float, default=0.0015,
                   help="Replicate-level change below this counts as rounding.")
    p.add_argument("--out", type=Path,
                   default=Path("paper/generated/structural-rebuild-diff.csv"))
    args = p.parse_args(argv)

    old = pd.read_csv(args.published).rename(columns={"structure_type": "state"})
    new = pd.read_csv(args.rebuilt).rename(columns={"structure_type": "state"})
    key = ["state", "genotype", "replicate"]
    metrics = {**DOR_METRICS, **RT_METRICS}

    # ---- replicate level ----------------------------------------------------
    m = old.merge(new, on=key, suffixes=("_old", "_new"), how="outer", indicator=True)
    print(f"replicate rows: published {len(old)}, rebuilt {len(new)}, "
          f"matched {int((m._merge == 'both').sum())}")
    if (m._merge != "both").any():
        print(m.loc[m._merge != "both", key + ["_merge"]].to_string(index=False))

    changed = []
    for label, col in metrics.items():
        if f"{col}_old" not in m.columns:
            continue
        a, b = m[f"{col}_old"], m[f"{col}_new"]
        if label in CIRCULAR:
            d = np.array([circular_delta(y, x) if np.isfinite(x) and np.isfinite(y)
                          else np.nan for x, y in zip(a, b)])
        else:
            d = (b - a).to_numpy(float)
        hit = np.isfinite(d) & (np.abs(d) > args.tol)
        for i in np.flatnonzero(hit):
            changed.append(dict(level="replicate", metric=label,
                                state=m.state.iat[i], genotype=m.genotype.iat[i],
                                replicate=int(m.replicate.iat[i]),
                                published=round(float(a.iat[i]), 3),
                                rebuilt=round(float(b.iat[i]), 3),
                                change=round(float(d[i]), 3)))
    rep = pd.DataFrame(changed)
    print(f"\nreplicate-level cells changed beyond {args.tol}: {len(rep)}")
    if len(rep):
        print(rep.groupby(["state", "genotype"]).size().to_string())

    # ---- Figure 5 level -----------------------------------------------------
    fig = []
    for label, col in metrics.items():
        if col not in old.columns or col not in new.columns:
            continue
        circ = label in CIRCULAR
        o = induced_fit(old, col, circ)
        n = induced_fit(new, col, circ)
        for geno in sorted(set(o.index) | set(n.index)):
            if geno not in o.index or geno not in n.index:
                continue
            ov, nv = o.loc[geno], n.loc[geno]
            delta = (circular_delta(nv.value, ov.value) if circ
                     else nv.value - ov.value)
            fig.append(dict(level="figure5", metric=label, genotype=geno,
                            published=round(float(ov.value), 2),
                            published_err=round(float(ov.err), 2),
                            rebuilt=round(float(nv.value), 2),
                            rebuilt_err=round(float(nv.err), 2),
                            change=round(float(delta), 2),
                            was_resolved=bool(ov.resolved),
                            now_resolved=bool(nv.resolved)))
    f5 = pd.DataFrame(fig)
    flipped = f5[f5.was_resolved != f5.now_resolved]
    print(f"\nFigure 5 cells: {len(f5)}; resolved before {int(f5.was_resolved.sum())}, "
          f"after {int(f5.now_resolved.sum())}; changed status {len(flipped)}")
    if len(flipped):
        print(flipped[["metric", "genotype", "published", "published_err",
                       "rebuilt", "rebuilt_err", "was_resolved",
                       "now_resolved"]].to_string(index=False))
    moved = f5.reindex(f5.change.abs().sort_values(ascending=False).index).head(12)
    print("\nlargest Figure 5 movements:")
    print(moved[["metric", "genotype", "published", "rebuilt", "change",
                 "was_resolved", "now_resolved"]].to_string(index=False))

    args.out.parent.mkdir(parents=True, exist_ok=True)
    pd.concat([rep, f5], ignore_index=True).to_csv(args.out, index=False)
    print(f"\nwrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
