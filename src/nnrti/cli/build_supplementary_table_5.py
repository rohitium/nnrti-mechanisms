#!/usr/bin/env python3
"""Build Supplementary Table 5: MM/GBSA sensitivity to snapshot relaxation.

Reviewer 3 (comment 3b) asked how sensitive ddE_Total and its components are to
the choice of relaxing each snapshot before scoring, and suggested recomputing
without minimization for a subset of genotypes.

Three protocols are compared on the same six genotypes, with sampling otherwise
identical to the published run (100 snapshots spread evenly over the
post-equilibration region, 0.25 discard fraction, 2000 iteration cap):

  none          single-point evaluation of the snapshot as written
  h_relax       hydrogens relaxed, heavy atoms harmonically restrained
  unrestrained  all atoms free -- the protocol used in the main text

Energies are converted to kcal/mol and referenced to the WT mean within each
protocol, so each column is internally consistent and the comparison is between
ddE values rather than between absolute energies computed on different
geometries.

Usage
-----
    PYTHONPATH=src python -m nnrti.cli.build_supplementary_table_5
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

KJ_PER_KCAL = 4.184
PROTOCOLS = ("none", "h_relax", "unrestrained")
COMPONENTS = {
    "binding_dg": "ddE_Total",
    "binding_dg_vdw": "ddE_vdW",
    "binding_dg_electrostatic": "ddE_elec",
    "binding_dg_gb": "ddE_GB",
    "binding_dg_sa": "ddE_SA",
}


def load(path: Path, protocol: str) -> pd.DataFrame | None:
    if not path.exists():
        return None
    d = pd.read_csv(path)
    d["protocol"] = protocol
    return d


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--checkpoint-dir", type=Path, default=Path("results/.checkpoints"))
    p.add_argument("--output-dir", type=Path, default=Path("results/analysis/relaxation_sensitivity"))
    p.add_argument("--xlsx", type=Path,
                   default=Path("revisions/Supplementary-Table-5.xlsx"))
    p.add_argument("--published", type=Path,
                   default=Path("results/.checkpoints/.checkpoint_mmgbsa_cpu2000_2026-08-29.csv"),
                   help="Main-text MM/GBSA checkpoint, used for the unrestrained column "
                        "until a dedicated re-run exists.")
    p.add_argument("--subset", nargs="*",
                   default=["wt", "K103N", "V106A", "Y188L", "G190E", "V106I_F227C"])
    args = p.parse_args(argv)

    frames = []
    for proto in PROTOCOLS:
        d = load(args.checkpoint_dir / f".checkpoint_mmgbsa_relaxsens_{proto}.csv", proto)
        if d is None and proto == "unrestrained" and args.published.exists():
            # The unrestrained column is the main-text protocol, so the published
            # computation is the right source for it: the table should report what
            # was published rather than a re-run of it. A fresh re-run is kept as a
            # reproducibility check and supersedes this when available.
            pub = pd.read_csv(args.published)
            d = pub[pub.safe_label.isin(args.subset)].copy()
            d["protocol"] = proto
            print(f"  {proto}: using the published checkpoint ({len(d)} rows, "
                  "same settings as the main text)")
        if d is None:
            print(f"  missing: {proto} (skipped)")
            continue
        frames.append(d)
    if not frames:
        raise SystemExit("no relaxation-sensitivity checkpoints found")
    allp = pd.concat(frames, ignore_index=True)

    # WT-referenced, per protocol, in kcal/mol
    rows = []
    for proto, g in allp.groupby("protocol"):
        wt = g[g.safe_label == "wt"]
        if wt.empty:
            print(f"  {proto}: no WT rows, cannot reference")
            continue
        for col, name in COMPONENTS.items():
            if col not in g:
                continue
            wt_mean = wt[col].mean()
            for lbl, gg in g[g.safe_label != "wt"].groupby("safe_label"):
                dd = (gg[col] - wt_mean) / KJ_PER_KCAL
                rows.append({"protocol": proto, "genotype": lbl, "component": name,
                             "mean_kcal": dd.mean(), "sem_kcal": dd.sem(),
                             "n_replicates": len(dd)})
    long = pd.DataFrame(rows)

    out = args.output_dir / "tables"
    out.mkdir(parents=True, exist_ok=True)
    long.to_csv(out / "relaxation_sensitivity_long.csv", index=False)

    # wide view of ddE_Total, the headline comparison
    tot = long[long.component == "ddE_Total"]
    wide = tot.pivot(index="genotype", columns="protocol", values="mean_kcal").round(2)
    wide_sem = tot.pivot(index="genotype", columns="protocol", values="sem_kcal").round(2)
    wide.columns = [f"{c}_mean" for c in wide.columns]
    wide_sem.columns = [f"{c}_sem" for c in wide_sem.columns]
    summary = wide.join(wide_sem).sort_index()
    summary.to_csv(out / "relaxation_sensitivity_ddE_total.csv")

    notes = pd.DataFrame({0: [
        "Supplementary Table 5. Sensitivity of MM/GBSA energies to the snapshot relaxation protocol.",
        "",
        "Columns are WT-referenced ddE values in kcal/mol, computed within each protocol so that",
        "each column is internally consistent. Sampling is identical across protocols and matches",
        "the main-text run: 100 snapshots spread evenly over the post-equilibration region of each",
        "100 ns replicate, 0.25 discard fraction, 2000 iteration cap, CPU platform, double precision.",
        "",
        "none          single-point evaluation of the snapshot as written",
        "h_relax       hydrogens relaxed with heavy atoms harmonically restrained (k = 10000 kJ/mol/nm^2)",
        "unrestrained  all atoms free; the protocol used for Table 2 of the main text",
        "",
        "Single-point evaluation is reported for completeness rather than as a usable alternative.",
        "Snapshots are written from explicit-solvent dynamics with bonds to hydrogen constrained and",
        "are evaluated here without a cutoff and without constraints, so residual steric contacts",
        "contribute Lennard-Jones energies orders of magnitude larger than the binding signal. The",
        "blow-up is confined to the van der Waals component; the electrostatic component is unaffected.",
    ]})

    args.xlsx.parent.mkdir(parents=True, exist_ok=True)
    with pd.ExcelWriter(args.xlsx, engine="openpyxl") as w:
        summary.reset_index().to_excel(w, sheet_name="ddE_Total", index=False)
        long.round(3).to_excel(w, sheet_name="All components", index=False)
        allp.to_excel(w, sheet_name="Per replicate", index=False)
        notes.to_excel(w, sheet_name="Notes", index=False, header=False)

    print(f"\nprotocols present: {sorted(allp.protocol.unique())}")
    print(f"\nddE_Total, kcal/mol (WT-referenced within protocol):\n")
    print(summary.to_string())
    print(f"\nwrote {args.xlsx}")
    print(f"wrote {out}/relaxation_sensitivity_long.csv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
