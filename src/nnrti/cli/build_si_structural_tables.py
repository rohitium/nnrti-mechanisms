#!/usr/bin/env python3
"""Build formula-based Supplementary Table 4 from full-precision replicate CSVs.

The workbook contains both-state inputs, summaries, within-state comparisons,
binding-associated interactions and a formula helper sheet. Replicate means
are the statistical unit; χ1 uses circular statistics. Excel or LibreOffice
must recalculate newly generated formula caches before display-only use.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats as sstats

from nnrti.cli.build_supplementary_table_4 import CATEGORY, PANEL_ORDER
from nnrti.cli.plot_structural_heatmap import (
    CIRCULAR,
    DOR_METRICS,
    RT_METRICS,
    SIGNIFICANCE_BY_PANEL,
    induced_fit_vs_wt,
    shift_vs_wt,
)

#: Decimals per metric, chosen so a column's smallest resolved difference is
#: still visible. Counts and volumes get one, distances and angles two.
_DECIMALS = {
    "chlorocyanophenyl_contacts_n": 1, "pyridinone_contacts_n": 1,
    "Tyr188_chlorophenyl_angle_deg": 1, "nnibp_volume_A3": 1,
    "chi1_Tyr181_deg": 1, "chi1_Tyr188_deg": 1,
}
_DEFAULT_DECIMALS = 2

#: State labels as they appear in the workbook, in the order they are listed.
_STATES = {"holo": "Bound (RT+DOR)", "apo": "Unliganded (RT only)"}


def _rank(genotypes: pd.Series) -> pd.Series:
    order = {m: i for i, m in enumerate(PANEL_ORDER)}
    return genotypes.map(lambda m: order.get(m, len(order)))


def _dec(col: str) -> int:
    return _DECIMALS.get(col, _DEFAULT_DECIMALS)


def _circ_mean(deg: np.ndarray) -> float:
    r = np.deg2rad(np.asarray(deg, float))
    r = r[np.isfinite(r)]
    if r.size == 0:
        return np.nan
    return float(np.rad2deg(np.arctan2(np.sin(r).mean(), np.cos(r).mean())))


def _circ_sem(deg: np.ndarray) -> float:
    r = np.deg2rad(np.asarray(deg, float))
    r = r[np.isfinite(r)]
    if r.size == 0:
        return np.nan
    R = float(np.hypot(np.sin(r).mean(), np.cos(r).mean()))
    sd = float(np.rad2deg(np.sqrt(-2.0 * np.log(max(R, 1e-12)))))
    return sd / np.sqrt(r.size)


def _fmt(mean: float, sem: float, decimals: int) -> str:
    """mean ± sem, widened when the error would otherwise print as zero.

    A cell reading "3.30 ± 0.00" claims an exact measurement. Where the error
    rounds away at the column's precision, show it to one significant figure
    instead, as the submitted Table 4 did.
    """
    if not np.isfinite(sem):
        return f"{mean:.{decimals}f}"
    d = decimals
    while d < decimals + 4 and sem > 0 and round(sem, d) == 0:
        d += 1
    return f"{mean:.{decimals}f} ± {sem:.{d}f}"


def _metric_columns() -> list[tuple[str, str]]:
    """(label, source column) in sheet order: drug-facing first, then RT's own."""
    return [(l, c) for l, c in DOR_METRICS.items()] + [(l, c) for l, c in RT_METRICS.items()]


def table4(per_rep: pd.DataFrame) -> dict[str, pd.DataFrame]:
    cols = _metric_columns()

    detail = pd.DataFrame({
        "Category": per_rep.genotype.map(CATEGORY),
        "Genotype": per_rep.genotype,
        "State": per_rep.structure_type.map(_STATES),
        "Replicate": per_rep.replicate,
    })
    for label, col in cols:
        detail[label] = pd.to_numeric(per_rep.get(col), errors="coerce")
    detail = (detail.assign(_r=_rank(detail.Genotype),
                            _s=detail.State.map({v: i for i, v in enumerate(_STATES.values())}))
              .sort_values(["_r", "_s", "Replicate"]).drop(columns=["_r", "_s"]))

    rows = []
    for (gt, st), grp in per_rep.groupby(["genotype", "structure_type"], sort=False):
        rec = {"Category": CATEGORY.get(gt, ""), "Genotype": gt, "State": _STATES[st],
               "Replicates": int(grp.replicate.nunique())}
        for label, col in cols:
            vals = pd.to_numeric(grp.get(col), errors="coerce").dropna()
            if vals.empty:
                rec[label] = ""
                continue
            if label in CIRCULAR:
                mean, sem = _circ_mean(vals.to_numpy()), _circ_sem(vals.to_numpy())
            else:
                mean = float(vals.mean())
                sem = float(vals.sem()) if vals.size > 1 else np.nan
            rec[label] = _fmt(mean, sem, _dec(col))
        rows.append(rec)
    summary = pd.DataFrame(rows)
    summary = (summary.assign(_r=_rank(summary.Genotype),
                              _s=summary.State.map({v: i for i, v in enumerate(_STATES.values())}))
               .sort_values(["_r", "_s"]).drop(columns=["_r", "_s"]))

    return {"Summary": summary.reset_index(drop=True),
            "Per-replicate": detail.reset_index(drop=True)}


def _pooled_df(sems, ns) -> float:
    """Sum n−1 over independent groups (df=4 or 8 for groups of three)."""
    del sems
    total = sum(n - 1 for n in ns if n >= 2)
    return float(total) if total > 0 else np.nan


def _group_stats(frame: pd.DataFrame, col: str, circular: bool) -> dict:
    """Per genotype: mean, standard error and n, on the circle where needed."""
    out = {}
    for gt, grp in frame.groupby("genotype"):
        v = pd.to_numeric(grp.get(col), errors="coerce").dropna()
        if v.empty:
            continue
        if circular:
            mean, sem = _circ_mean(v.to_numpy()), _circ_sem(v.to_numpy())
        else:
            mean = float(v.mean())
            sem = float(v.sem()) if v.size > 1 else np.nan
        out[gt] = (mean, sem, int(v.size))
    return out


def _long(val, err, tcrit: float, state: str, components: dict) -> pd.DataFrame:
    """One row per genotype and metric: the figure's cell, and the test behind it.

    `components` carries the per-group standard errors and counts that the cell
    was built from, which is what the degrees of freedom need. Deriving them
    here rather than assuming 4 and 8 means a genotype with a missing replicate
    gets its own df instead of the panel's nominal one.
    """
    rows = []
    for gt in val.index:
        for metric in val.columns:
            v, e = val.loc[gt, metric], err.loc[gt, metric]
            if not np.isfinite(v):
                continue
            sems, ns = components.get((gt, metric), ((), ()))
            df = _pooled_df(sems, ns) if ns else np.nan
            half = float(sstats.t.ppf(0.975, df)) * e if np.isfinite(e) and np.isfinite(df) else np.nan
            t = float(v) / float(e) if np.isfinite(e) and e > 0 else np.nan
            p = (float(2 * sstats.t.sf(abs(t), df))
                 if np.isfinite(t) and np.isfinite(df) else np.nan)
            rows.append({
                "Category": CATEGORY.get(gt, ""), "Genotype": gt, "State": state,
                "Metric": metric,
                "Shift": round(float(v), 3),
                "SEM": round(float(e), 3) if np.isfinite(e) else "",
                "t": round(t, 2) if np.isfinite(t) else "",
                "df": int(df) if np.isfinite(df) else "",
                "p": (float(f"{p:.3g}") if np.isfinite(p) else ""),
                "CI95 low": round(float(v - half), 3) if np.isfinite(half) else "",
                "CI95 high": round(float(v + half), 3) if np.isfinite(half) else "",
                "95% CI excludes zero": ("yes" if np.isfinite(half) and abs(v) > half else "no"),
            })
    out = pd.DataFrame(rows)
    return out.assign(_r=_rank(out.Genotype)).sort_values(["_r", "Metric"]).drop(columns="_r")


def statistics_sheets(per_rep: pd.DataFrame) -> dict[str, pd.DataFrame]:
    """The two comparisons the figures draw, each with the test behind it.

    One sheet per panel, because the two use different degrees of freedom and
    therefore different 95% criteria, and putting them in one sheet invites a
    reader to apply the wrong one.
    """
    frames = []
    for state, label in _STATES.items():
        metrics = {**DOR_METRICS, **RT_METRICS} if state == "holo" else RT_METRICS
        val, err = shift_vs_wt(per_rep, state, metrics)
        if val.empty:
            continue
        # Two group means are subtracted here: the genotype's and WT's.
        comp = {}
        frame = per_rep[per_rep.structure_type == state]
        for metric, col in metrics.items():
            if metric not in val.columns:
                continue
            g = _group_stats(frame, col, metric in CIRCULAR)
            if "WT" not in g:
                continue
            for gt, (_m, sem, n) in g.items():
                comp[(gt, metric)] = ((sem, g["WT"][1]), (n, g["WT"][2]))
        frames.append(_long(val, err, SIGNIFICANCE_BY_PANEL["A"], label, comp))
    mutant = pd.concat(frames, ignore_index=True)

    val, err = induced_fit_vs_wt(per_rep, RT_METRICS)
    # Four group means are subtracted: holo and apo, for the genotype and WT.
    comp = {}
    holo = per_rep[per_rep.structure_type == "holo"]
    apo = per_rep[per_rep.structure_type == "apo"]
    for metric, col in RT_METRICS.items():
        if metric not in val.columns:
            continue
        h, a = _group_stats(holo, col, metric in CIRCULAR), _group_stats(apo, col, metric in CIRCULAR)
        if "WT" not in h or "WT" not in a:
            continue
        for gt in set(h) & set(a):
            comp[(gt, metric)] = ((h[gt][1], a[gt][1], h["WT"][1], a["WT"][1]),
                                  (h[gt][2], a[gt][2], h["WT"][2], a["WT"][2]))
    induced = (_long(val, err, SIGNIFICANCE_BY_PANEL["B"], "Bound minus unliganded", comp)
               .drop(columns="State").reset_index(drop=True))

    # The 95% column comes from the panel's nominal criterion and p comes from
    # the row's own degrees of freedom. They agree unless a genotype is short a
    # replicate, so a disagreement is a signal, not noise: report it rather
    # than let the two columns quietly contradict each other.
    for name, sheet in (("Mutant vs WT", mutant), ("Induced fit vs WT", induced)):
        num = pd.to_numeric(sheet["p"], errors="coerce")
        clash = ((num < 0.05) != (sheet["95% CI excludes zero"] == "yes")) & num.notna()
        if clash.any():
            print(f"  {name}: {clash.sum()} row(s) where p and the panel criterion differ "
                  f"(a genotype with fewer than three usable replicates, so its df is "
                  f"lower than the panel's): "
                  + ", ".join(f"{r.Genotype}/{r.Metric}" for r in sheet[clash].itertuples()))
    return {"Mutant vs WT": mutant, "Induced fit vs WT": induced}


def _write(path: Path, sheets: dict[str, pd.DataFrame]) -> None:
    with pd.ExcelWriter(path, engine="openpyxl") as xl:
        for name, df in sheets.items():
            df.to_excel(xl, sheet_name=name, index=False)
            ws = xl.sheets[name]
            ws.freeze_panes = "A2"
            for i, col in enumerate(df.columns, start=1):
                width = max(len(str(col)), *(len(str(v)) for v in df[col].head(80))) + 2
                ws.column_dimensions[ws.cell(row=1, column=i).column_letter].width = min(width, 72)
    print(f"wrote {path}  sheets: {list(sheets)}")


def main(argv: list[str] | None = None) -> int:
    root = Path(__file__).resolve().parents[3]
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--per-replicate", type=Path,
                   default=root / "results/analysis/structural_metrics/tables/"
                                  "structural_metrics_per_replicate.csv")
    p.add_argument("--out-dir", type=Path, default=root / "paper/generated/tables")
    a = p.parse_args(argv)

    per_rep = pd.read_csv(a.per_replicate)
    a.out_dir.mkdir(parents=True, exist_ok=True)
    sheets = table4(per_rep)
    stats = statistics_sheets(per_rep)
    output = a.out_dir / "Supplementary-Table-4.xlsx"
    _write(output, {**sheets, **stats})
    from nnrti.analysis.structural_workbook import add_formulas
    add_formulas(output)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
