#!/usr/bin/env python3
"""Structural metrics for both legs: the data behind Supplementary Table 4.

Reviewer 1 (comment 4) asked for protein-level observables the manuscript did
not show -- RMSF, secondary-structure changes and pocket conformational
distributions -- to separate direct changes in DOR contacts from broader
remodelling of the pocket, and (comment 3) for the same comparison in apo and
holo RT.

Every metric that can be measured without a ligand is reported for both legs,
tagged by ``structure_type``. The DOR-contact metrics are necessarily holo-only
and are left empty for apo rather than dropped, so a reader can see at a glance
which columns the question even applies to.

Usage
-----
    PYTHONPATH=src python -m nnrti.cli.build_structural_metrics_table
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

RESID_OFFSET = -3

#: Frames with a pocket Ca radius of gyration above this are excluded, matching
#: compare_apo_ensembles.py and compute_protein_level_changes.py so that every
#: apo number in the paper rests on the same frames.
POCKET_RG_MAX_A = 10.0

NNIBP_P66 = (100, 101, 103, 106, 107, 108, 179, 181, 188, 189, 190, 227, 229, 234, 318)

#: Column-safe region names, matching compute_protein_level_changes.REGIONS.
REGIONS = {
    "nnibp_lining": NNIBP_P66,
    "primer_grip": tuple(range(225, 236)),
    "beta6": tuple(range(100, 111)),
    "beta9_beta10": tuple(range(179, 191)),
    "ymdd": (183, 184, 185, 186),
    "fingers": tuple(range(1, 86)) + tuple(range(118, 156)),
    "palm": tuple(range(86, 118)) + tuple(range(156, 237)),
    "thumb": tuple(range(237, 319)),
    "connection": tuple(range(319, 427)),
}

#: DOR-contact metrics. Holo only: there is no ligand in the apo leg.
CONTACT_METRICS = {
    "res103_bb_to_triazN": "res103_CO_to_triazolinone_N_A",
    "res103_sc_to_dor_polar": "res103_sidechain_to_DOR_A",
    "chl_ring_burial": "chlorocyanophenyl_contacts_n",
    "y188_interplanar_deg": "Tyr188_chlorophenyl_angle_deg",
    "res179_to_dor_mindist": "res179_to_DOR_A",
    "dor_to_S105_mindist": "Ser105_to_DOR_A",
    "dor_to_res106_mindist": "res106_to_DOR_A",
}


def _display(label: str) -> str:
    return "WT" if label.lower() == "wt" else label.upper().replace("_", "+")


def _regional(per_residue: Path, value_cols: dict[str, str]) -> pd.DataFrame:
    """Mean of each per-residue quantity over each region, per replicate (p66)."""
    d = pd.read_csv(per_residue)
    d = d[d.chain == 0]
    out = []
    for (st, lbl, rep), g in d.groupby(["structure_type", "safe_label", "replicate"]):
        row = {"structure_type": st, "genotype": _display(lbl), "replicate": rep}
        for region, resids in REGIONS.items():
            sel = g[g.canonical.isin(resids)]
            for src, suffix in value_cols.items():
                row[f"{suffix}_{region}"] = (round(float(sel[src].mean()), 3)
                                             if len(sel) else np.nan)
        out.append(row)
    return pd.DataFrame(out)


def build(root: Path) -> pd.DataFrame:
    # ---- protein-level, both legs -------------------------------------------
    frames = []
    for state, base in (("holo", "protein_level"), ("apo", "protein_level_apo")):
        t = root / f"results/analysis/{base}/tables"
        rmsf = _regional(t / "ca_rmsf_per_residue.csv", {"rmsf": "rmsf_A"})
        ss = _regional(t / "secondary_structure_per_residue.csv",
                       {"frac_helix": "helix_frac", "frac_strand": "strand_frac"})
        frames.append(rmsf.merge(ss, on=["structure_type", "genotype", "replicate"],
                                 how="outer"))
    out = pd.concat(frames, ignore_index=True)

    # ---- chi1 of the two pocket tyrosines, both legs ------------------------
    chi = pd.read_csv(root / "results/analysis/chi1_apo_holo/tables/chi1_per_replicate.csv")
    chi["genotype"] = chi.genotype.map(_display)
    chi = chi.rename(columns={"state": "structure_type"})
    wide = chi.pivot_table(index=["structure_type", "genotype", "replicate"],
                           columns="residue", values="chi1_mean").reset_index()
    wide = wide.rename(columns={"Tyr181": "chi1_Tyr181_deg", "Tyr188": "chi1_Tyr188_deg"})
    out = out.merge(wide, on=["structure_type", "genotype", "replicate"], how="left")

    # ---- NNIBP volume, both legs --------------------------------------------
    holo_v = pd.read_csv(root / "results/analysis/pocket_volume/tables/pocket_volume_per_rep.csv")
    holo_v = pd.DataFrame({"structure_type": "holo", "genotype": holo_v.mutation,
                           "replicate": holo_v.replicate,
                           "nnibp_volume_A3": holo_v.pocket_volume_mean.round(1)})
    apo_f = pd.read_csv(root / "results/analysis/apo_ensembles/tables/apo_frame_metrics.csv")
    # Drop frames whose pocket has come apart, at the threshold
    # compare_apo_ensembles.py and compute_protein_level_changes.py both use.
    # Without it, apo g190e rep_02 averages 1871 A^3 (max 4172) instead of 103,
    # because a disintegrated pocket has an enormous convex hull -- which
    # reported G190E's apo pocket as +575 A^3 against WT.
    # The same three filters compare_apo_ensembles and the opening cost apply,
    # in the same order. Using a different subset here is how three values for
    # the same apo pocket volume ended up in circulation.
    before = len(apo_f)
    if "pocket_intact" in apo_f.columns:
        apo_f = apo_f[apo_f.pocket_intact]
    else:
        apo_f = apo_f[apo_f.pocket_rg <= POCKET_RG_MAX_A]
    apo_f = apo_f[apo_f.window == "matched"]
    total = apo_f.groupby(["safe_label", "replicate"]).total_ns.first()
    window_ns = apo_f.time_ns.max()
    keep = {k for k, v in total.items() if v >= window_ns * 0.98}
    dropped_reps = sorted(set(total.index) - keep)
    short_apo = {(_display(l), r) for l, r in dropped_reps}
    apo_f = apo_f[[(l, r) in keep for l, r in zip(apo_f.safe_label, apo_f.replicate)]]
    if before != len(apo_f):
        print(f"  apo: kept {len(apo_f)} of {before} frames")
        if dropped_reps:
            print("  apo: excluded replicates shorter than the matched window: "
                  + ", ".join(f"{l} rep{r}" for l, r in dropped_reps))
    apo_v = (apo_f.groupby(["safe_label", "replicate"], as_index=False)
             .pocket_volume.mean())
    apo_v = pd.DataFrame({"structure_type": "apo",
                          "genotype": apo_v.safe_label.map(_display),
                          "replicate": apo_v.replicate,
                          "nnibp_volume_A3": apo_v.pocket_volume.round(1)})
    out = out.merge(pd.concat([holo_v, apo_v], ignore_index=True),
                    on=["structure_type", "genotype", "replicate"], how="left")

    # ---- DOR contacts: holo only --------------------------------------------
    # mechanism_coordinates.csv is PER FRAME (139-540 rows per replicate), so it
    # must be averaged to one row per replicate before merging. Merging it raw
    # multiplies every protein-level row by the frame count: the first version
    # of this produced 14,769 holo rows instead of 60.
    coords = pd.read_csv(root / "results/analysis/mechanisms/mechanism_coordinates.csv")
    moiety = pd.read_csv(root / "results/analysis/mechanisms/dor_moiety_contacts_per_replicate.csv")
    c = (coords.groupby(["mutation", "replicate"], as_index=False)[list(CONTACT_METRICS)]
         .mean())
    c = c.merge(moiety[["mutation", "replicate", "pyridinone"]],
                on=["mutation", "replicate"], how="left")
    c = c.rename(columns={**CONTACT_METRICS, "pyridinone": "pyridinone_contacts_n",
                          "mutation": "genotype"})
    c.insert(0, "structure_type", "holo")
    out = out.merge(c, on=["structure_type", "genotype", "replicate"], how="left")

    # A replicate too short for the matched window cannot take part in ANY
    # matched apo comparison, not just the volume one. Leaving its RMSF and
    # secondary structure in would compare a 34 ns run against 100 ns ones.
    if short_apo:
        before_rows = len(out)
        mask = pd.Series(
            [(st == "apo" and (g, r) in short_apo)
             for st, g, r in zip(out.structure_type, out.genotype, out.replicate)],
            index=out.index)
        out = out[~mask]
        print(f"  apo: dropped {before_rows - len(out)} replicate row(s) too "
              f"short for the matched window: "
              + ", ".join(f"{g} rep{r}" for g, r in sorted(short_apo)))

    lead = ["structure_type", "genotype", "replicate"]
    contact_cols = [v for v in CONTACT_METRICS.values()] + ["pyridinone_contacts_n"]
    pocket_cols = ["nnibp_volume_A3", "chi1_Tyr181_deg", "chi1_Tyr188_deg"]
    rest = [c for c in out.columns if c not in lead + contact_cols + pocket_cols]
    return out[lead + contact_cols + pocket_cols + rest].sort_values(
        ["structure_type", "genotype", "replicate"]).reset_index(drop=True)


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--root", type=Path, default=Path("."))
    p.add_argument("--output-dir", type=Path,
                   default=Path("results/analysis/structural_metrics/tables"))
    args = p.parse_args(argv)

    per_rep = build(args.root)

    num = per_rep.select_dtypes("number").columns.drop("replicate")
    summary = (per_rep.groupby(["structure_type", "genotype"], as_index=False)
               .agg(n_replicates=("replicate", "nunique"),
                    **{c: (c, "mean") for c in num}))
    for c in num:
        summary[c] = summary[c].round(3)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    per_rep.to_csv(args.output_dir / "structural_metrics_per_replicate.csv", index=False)
    summary.to_csv(args.output_dir / "structural_metrics_summary.csv", index=False)

    # One row per leg, genotype and replicate. A merge against a per-frame table
    # inflates this silently, so it is checked rather than eyeballed.
    dup = per_rep.duplicated(["structure_type", "genotype", "replicate"]).sum()
    if dup:
        raise AssertionError(
            f"{dup} duplicated (structure_type, genotype, replicate) rows: a source "
            "table is per-frame and was not aggregated before merging")

    print(f"per-replicate: {per_rep.shape[0]} rows x {per_rep.shape[1]} cols")
    print(f"summary      : {summary.shape[0]} rows x {summary.shape[1]} cols")
    print(f"  structure types: {dict(per_rep.groupby('structure_type').size())}")
    holo_only = [c for c in CONTACT_METRICS.values()] + ["pyridinone_contacts_n"]
    apo = per_rep[per_rep.structure_type == "apo"]
    filled = {c: int(apo[c].notna().sum()) for c in holo_only}
    print(f"  DOR-contact columns filled in apo rows (should be 0): "
          f"{max(filled.values())}")
    print(f"\nwrote {args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
