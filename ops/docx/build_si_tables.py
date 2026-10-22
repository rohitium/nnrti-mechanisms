#!/usr/bin/env python3
"""Emit the Supporting Information tables as CSV, from committed analysis output.

Each table is re-derived here rather than typed, so a rerun of the underlying
analysis propagates into the document on the next insert instead of silently
disagreeing with it. Written to revisions/si_tables/.

Usage
-----
    PYTHONPATH=src python ops/docx/build_si_tables.py
"""

from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "revisions/si_tables"

#: apo safe_label -> the genotype name used in the main-text tables
DISPLAY = {
    "wt": "WT", "a98g_f227c": "A98G+F227C", "f227c": "F227C", "g190e": "G190E",
    "k103n_m230l": "K103N+M230L", "v106a": "V106A", "v106a_p225h": "V106A+P225H",
    "v106i_f227c": "V106I+F227C",
}


def charge_model() -> pd.DataFrame:
    d = pd.read_csv(ROOT / "results/analysis/relaxation_sensitivity/tables/"
                    "charge_model_comparison.csv")
    order = ["ddE_Total", "ddE_vdW", "ddE_elec", "ddE_GB"]
    d["component"] = pd.Categorical(d.component, order, ordered=True)
    d = d.sort_values(["genotype", "component"])
    rows = [["Genotype", "Component", "Gasteiger", "AM1-BCC", "Difference"]]
    for _, r in d.iterrows():
        rows.append([r.genotype, r.component.replace("ddE_", "ΔΔE "),
                     f"{r.mean_gast:.2f} ± {r.sem_gast:.2f}",
                     f"{r.mean_am1bcc:.2f} ± {r.sem_am1bcc:.2f}",
                     f"{r.difference:+.2f}"])
    return pd.DataFrame(rows[1:], columns=rows[0])


def apo_pocket() -> pd.DataFrame:
    a = pd.read_csv(ROOT / "results/analysis/apo_ensembles/tables/"
                    "apo_matched_window_summary.csv")
    t3 = pd.read_csv(ROOT / "paper/tables/Table-3-structural.csv")
    holo = {}
    for _, r in t3.iterrows():
        v = str(r["V(NNIBP)"])
        holo[r.Genotype] = float(v.split("±")[0]) if "±" in v else float(v)
    rows = [["Genotype", "Replicates", "Apo V(NNIBP) (Å³)",
             "Holo V(NNIBP) (Å³)", "Induced-fit gap (Å³)", "ΔApo vs WT (Å³)"]]
    for _, r in a.sort_values("pocket_volume").iterrows():
        name = DISPLAY.get(r.safe_label, r.safe_label)
        h = holo.get(name)
        rows.append([
            name, int(r.n_rep), f"{r.pocket_volume:.1f} ± {r.pocket_volume_sd:.1f}",
            f"{h:.1f}" if h is not None else "n/a",
            f"{h - r.pocket_volume:.0f}" if h is not None else "n/a",
            f"{r.dVolume_vs_WT:+.1f}"])
    return pd.DataFrame(rows[1:], columns=rows[0])


def relaxation() -> pd.DataFrame:
    d = pd.read_csv(ROOT / "results/analysis/relaxation_sensitivity/tables/"
                    "relaxation_sensitivity_long.csv")
    order = ["ddE_Total", "ddE_vdW", "ddE_elec", "ddE_GB"]
    d["component"] = pd.Categorical(d.component, order, ordered=True)
    rows = [["Genotype", "Component", "Single-point", "H-only relaxation",
             "Full minimization"]]
    for (g, c), sub in d.groupby(["genotype", "component"], observed=True):
        get = lambda p: sub[sub.protocol == p]
        def fmt(p, sci=False):
            s = get(p)
            if s.empty:
                return "n/a"
            m, e = float(s.mean_kcal.iloc[0]), float(s.sem_kcal.iloc[0])
            return f"{m:.2e} ± {e:.1e}" if abs(m) >= 1e4 else f"{m:.2f} ± {e:.2f}"
        rows.append([g.replace("_", "+"), str(c).replace("ddE_", "ΔΔE "),
                     fmt("none"), fmt("h_relax"), fmt("unrestrained")])
    return pd.DataFrame(rows[1:], columns=rows[0])


def relaxation_panel() -> pd.DataFrame:
    """Replacement Table 5: both relaxed protocols across the whole panel."""
    d = pd.read_csv(ROOT / "results/analysis/relaxation_sensitivity/tables/"
                    "relaxation_panel_long.csv")
    order = ["ddE_Total", "ddE_vdW", "ddE_elec", "ddE_GB"]
    d["component"] = pd.Categorical(d.component, order, ordered=True)
    d = d.sort_values(["genotype", "component"])
    rows = [["Genotype", "Component", "H-only relaxation", "Full minimization"]]
    for _, r in d.iterrows():
        rows.append([r.genotype, str(r.component).replace("ddE_", "\u0394\u0394E "),
                     f"{r.h_relax_mean:.2f} \u00b1 {r.h_relax_sem:.2f}",
                     f"{r.full_min_mean:.2f} \u00b1 {r.full_min_sem:.2f}"])
    return pd.DataFrame(rows[1:], columns=rows[0])


def relaxation_correlations() -> pd.DataFrame:
    """Agreement between protocols, and each protocol against phenotype."""
    a = pd.read_csv(ROOT / "results/analysis/relaxation_sensitivity/tables/"
                    "relaxation_protocol_agreement.csv")
    ph = pd.read_csv(ROOT / "results/analysis/relaxation_sensitivity/tables/"
                     "relaxation_phenotype_correlation.csv")
    rows = [["Component", "Protocols agree (Spearman)", "Sign agreement",
             "vs phenotype (H-only)", "vs phenotype (full min)"]]
    for _, r in a.iterrows():
        h = ph[(ph.component == r.component) & (ph.protocol == "h_relax")].iloc[0]
        f = ph[(ph.component == r.component) & (ph.protocol == "full_min")].iloc[0]
        rows.append([str(r.component).replace("ddE_", "\u0394\u0394E "),
                     f"{r.protocol_spearman:.2f} (p = {r.protocol_p:.3f})",
                     r.sign_agreement,
                     f"{h.spearman:+.2f}", f"{f.spearman:+.2f}"])
    return pd.DataFrame(rows[1:], columns=rows[0])


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    for name, fn in (("charge_model", charge_model),
                     ("apo_pocket", apo_pocket),
                     ("relaxation", relaxation),
                     ("relaxation_panel", relaxation_panel),
                     ("relaxation_correlations", relaxation_correlations)):
        df = fn()
        p = OUT / f"SI-{name}.csv"
        df.to_csv(p, index=False)
        print(f"wrote {p.relative_to(ROOT)}  ({len(df)} rows x {len(df.columns)} cols)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
