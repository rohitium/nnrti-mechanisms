#!/usr/bin/env python3
"""Intra-protein stand-ins for the protein-ligand metrics of Table 3.

Seven of the nine Table 3 columns measure a protein-ligand contact. They
therefore have no counterpart in the unliganded enzyme, which is why the
apo/holo comparison the reviewer asks for cannot be run on them directly:
there is no X(mutant, apo) to subtract.

This looks for a replacement that does exist in both states. A contact such as
d(Lys103 C=O - DOR triazolinone N) lies along a direction across the pocket.
Extend that direction past the ligand and it strikes the far wall. The distance
from the anchor atom to that far wall is the width of the pocket in the
direction of the contact -- an intra-protein distance, measurable with or
without a ligand, and expected to track the contact it was built from.

So, per Table 3 column:

    A = anchor atom group on the protein (the residue the column names)
    L = the ligand group it contacts
    u = unit vector A -> L
    B = protein heavy atoms lying beyond L along u, within a narrow cylinder

and the candidate surrogate is d(A, B).

Candidates are ranked by how well they track the Table 3 column *across
genotypes*, which is the axis the final test uses. A within-trajectory
correlation is thermal coupling and says nothing about whether a mutation moves
both together, so it is reported separately and never used to choose.

Selection never looks at ddG_bind. That keeps the later test honest: the
surrogate is chosen on its relation to the Table 3 column alone, so the one
test per column against ddG_bind is not a search over thousands of pairs.

A column for which no candidate tracks it is an informative failure, not a
missing result: it means the column reports where the ligand sat rather than
the shape of the pocket, so it has no apo analogue even in principle.

Usage
-----
    PYTHONPATH=src python -m nnrti.cli.find_pocket_surrogates --stage axes
    PYTHONPATH=src python -m nnrti.cli.find_pocket_surrogates --stage scan
    PYTHONPATH=src python -m nnrti.cli.find_pocket_surrogates --stage rank
"""
from __future__ import annotations

import argparse
import logging
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

from nnrti.cli.compare_apo_ensembles import frame_times_ns, true_steps_completed
from nnrti.cli.compute_mechanism_coordinates import (
    LIG_RESNAME,
    OFFSET,
    _ligand_rings,
    _mindist,
    _replicate_inputs,
    _res_atoms,
)

REPO_ROOT = Path(__file__).resolve().parents[3]

# Each entry: Table 3 column, the per-frame name it has in mechanism_coordinates,
# the anchor residue and atom names (None = whole side chain), and the ligand group.
AXES = [
    ("d(103 C=O...N4x)", "res103_bb_to_triazN", 103, {"N", "O"}, "triaz_N"),
    ("d(103 sc...DOR)", "res103_sc_to_dor", 103, None, "all_heavy"),
    ("d(V179...DOR)", "res179_to_dor_mindist", 179, None, "all_heavy"),
    ("d(S105...DOR)", "dor_to_S105_mindist", 105, None, "all_heavy"),
    ("d(106...DOR)", "dor_to_res106_mindist", 106, None, "all_heavy"),
]

# cylinder around the contact direction, in Angstrom
PERP_MAX = 4.0        # how far off the axis an atom may sit
REACH_MIN = 1.0       # how far beyond the ligand atom it must lie
REACH_MAX = 11.0      # and no further

# High fold resistance that ddG_bind does not reproduce -- the genotypes this
# whole comparison exists to explain.
MISMATCH = {"K103N+M230L", "A98G+F227C", "G190E"}

# Partner atoms are backbone only. The far wall of the pocket is lined by
# Phe227, Leu234, Tyr318 and Leu100, every one of which is mutated somewhere in
# this panel. A side-chain partner atom would then move because the side chain
# was replaced -- F227C deletes most of the ring -- and the surrogate would
# report chemistry rather than pocket shape. Backbone atoms exist in every
# mutant and mean the same thing in all of them.
BACKBONE_ATOMS = {"N", "CA", "C", "O"}


APO_ROOT = Path("results/md_runs/apo")


def _mic(delta: np.ndarray, box: np.ndarray) -> np.ndarray:
    """Minimum-image displacement, for a rectangular box.

    make_molecules_whole is not safe for these runs. It unwraps along the bond
    graph, so a single bond that happens to span more than half the box
    propagates its error through the rest of the chain: on apo G190E rep_02 it
    reports d(Lys103, Pro226) as 47 A where the true distance is 14 A, which is
    what made that replicate look unfolded. Every distance used here is far
    below half the box edge (about 14-20 A against 64 A), so the minimum image
    is unambiguous and needs no unwrapping at all.
    """
    return delta - box * np.round(delta / box)


def _group_centroid(xyz: np.ndarray, idx: list[int], box: np.ndarray) -> np.ndarray:
    """Centroid of an atom group, imaged against its own first atom."""
    sub = xyz[:, idx, :]
    ref = sub[:, :1, :]
    return ref[:, 0, :] + _mic(sub - ref, box[:, None, :]).mean(axis=1)


def _run_ns(rep_dir: Path) -> float:
    """Nanoseconds actually on disk for a replicate, from its state CSV."""
    js = next((q for q in sorted(rep_dir.glob("*.json"))), None)
    claimed = 0
    if js is not None:
        import json as _json
        try:
            claimed = int(_json.loads(js.read_text())
                          .get("md_production_steps_completed") or 0)
        except Exception:
            claimed = 0
    return true_steps_completed(rep_dir, claimed) * 2.0 / 1_000_000.0


def _window_mask(n_frames: int, total_ns: float, last_ns: float) -> np.ndarray:
    """Frames inside the final `last_ns` of a run.

    The DCD schedule changes stride when a run is extended, so the frames are
    not evenly spaced in time and a plain mean over them over-weights the
    densely sampled head. Selecting a matched trailing window in nanoseconds
    makes apo and holo means comparable across runs of different length.
    """
    t = frame_times_ns(n_frames, total_ns)
    return t >= (t[-1] - last_ns)

# holo manifest labels -> apo directory names
def _apo_dirname(mutation: str) -> str:
    return mutation.lower().replace("+", "_")


def _apo_replicates(mutation: str) -> list[tuple[int, Path, Path]]:
    """(replicate, topology, trajectory) for every apo replicate of a genotype.

    The apo runs are not in the holo manifest, so they are discovered on disk.
    Preference goes to the periodicity-corrected trajectory where one exists;
    make_molecules_whole is applied either way, so both are safe.
    """
    out = []
    for d in sorted((REPO_ROOT / APO_ROOT).glob(f"{_apo_dirname(mutation)}/rep_*")):
        # Always the raw trajectory. 19 of 45 apo *_pbcfix.dcd files disagree in
        # frame count with the run they claim to correct -- some are longer than
        # their source -- so they are leftovers from superseded runs. Periodicity
        # is handled here by minimum image, which makes them unnecessary anyway.
        raw = [q for q in d.glob("*_analysis.dcd") if "pbcfix" not in q.name]
        topo = [q for q in d.glob("*_analysis_topology.pdb") if "from_md" not in q.name]
        src = sorted(raw)
        if not src or not topo:
            continue
        out.append((int(d.name.split("_")[-1]), sorted(topo)[-1], src[-1]))
    return out


def _anchor(top, canonical: int, names):
    return _res_atoms(top, canonical, sidechain_only=names is None, names=names)


def _candidates(traj, anchor_idx, lig_idx, canonical: int,
                backbone_only: bool = True) -> pd.DataFrame:
    """Protein atoms lying beyond the ligand group along the contact axis."""
    top = traj.topology
    exclude = {canonical + OFFSET + d for d in (-2, -1, 0, 1, 2)}
    pool = [
        a.index for a in top.atoms
        if a.residue.name != LIG_RESNAME
        and a.element.symbol != "H"
        and a.residue.resSeq not in exclude
        and (a.name in BACKBONE_ATOMS if backbone_only else True)
    ]
    xyz = traj.xyz * 10.0                      # (frames, atoms, 3)
    A = xyz[:, anchor_idx, :].mean(axis=1)     # (frames, 3)
    L = xyz[:, lig_idx, :].mean(axis=1)
    axis = L - A
    span = np.linalg.norm(axis, axis=1)        # (frames,)
    u = axis / span[:, None]

    rel = xyz[:, pool, :] - A[:, None, :]      # (frames, pool, 3)
    proj = np.einsum("fpc,fc->fp", rel, u)
    perp = np.linalg.norm(rel - proj[:, :, None] * u[:, None, :], axis=2)
    reach = proj - span[:, None]               # distance beyond the ligand group

    keep = (
        (perp.mean(axis=0) <= PERP_MAX)
        & (reach.mean(axis=0) >= REACH_MIN)
        & (reach.mean(axis=0) <= REACH_MAX)
    )
    rows = []
    for j, idx in enumerate(pool):
        if not keep[j]:
            continue
        atom = top.atom(idx)
        rows.append(dict(
            canonical_anchor=canonical,
            partner_resid=atom.residue.resSeq,
            partner_canonical=atom.residue.resSeq - OFFSET,
            partner_resname=atom.residue.name,
            partner_atom=atom.name,
            mean_perp_angstrom=round(float(perp[:, j].mean()), 2),
            mean_reach_angstrom=round(float(reach[:, j].mean()), 2),
        ))
    return pd.DataFrame(rows)


def _anchor_variants(top, canonical: int, names):
    """The two anchor groups a surrogate may be measured from.

    "as_table3" uses the same atoms the Table 3 column does, so the surrogate is
    directly comparable to it. "ca" uses the anchor residue's CA alone. Residues
    103 and 106 are mutated in this panel, so their side chains differ between
    genotypes; a CA anchor is identity-independent and reports where the residue
    sits rather than what it is made of.
    """
    out = {}
    group = _anchor(top, canonical, names)
    if group:
        out["as_table3"] = group
    ca = _res_atoms(top, canonical, names={"CA"})
    if ca:
        out["ca"] = ca
    return out


def stage_scan(args) -> int:
    import mdtraj as md

    cand = pd.read_csv(args.candidates)
    mf = pd.read_csv(args.manifest)
    if args.mutations:
        mf = mf[mf["mutation"].isin(args.mutations)]
    specs = {c: dict(canonical=a, names=n) for c, co, a, n, lg in
             [(c, co, a, n, lg) for c, co, a, n, lg in AXES]}

    jobs = []
    for _, r in mf.iterrows():
        mutation = str(r["mutation"])
        if args.state in ("holo", "both"):
            try:
                topo, dcd, _ = _replicate_inputs(r, REPO_ROOT)
                jobs.append(("holo", mutation, int(r["replicate"]), topo, dcd))
            except Exception as exc:
                logging.error("holo %s rep%s: %s", mutation, r["replicate"], exc)
    if args.state in ("apo", "both"):
        for mutation in dict.fromkeys(str(m) for m in mf["mutation"]):
            for rep, topo, dcd in _apo_replicates(mutation):
                jobs.append(("apo", mutation, rep, topo, dcd))

    rows = []
    for state, mutation, rep, topo, dcd in jobs:
        try:
            traj = md.load(str(dcd), top=str(topo), stride=args.stride)
            top = traj.topology
            xyz = traj.xyz * 10.0
            if traj.unitcell_lengths is None:
                raise ValueError("no box information; minimum image is undefined")
            box = traj.unitcell_lengths * 10.0

            per_frame = {}
            for column, grp in cand.groupby("table3_column", sort=False):
                spec = specs[column]
                anchors = _anchor_variants(top, spec["canonical"], spec["names"])
                for _, c in grp.iterrows():
                    partner = [
                        a.index for a in top.atoms
                        if a.residue.resSeq == c.partner_resid and a.name == c.partner_atom
                    ]
                    if not partner:
                        continue
                    for tag, anchor_idx in anchors.items():
                        centre = _group_centroid(xyz, list(anchor_idx), box)
                        d = np.linalg.norm(
                            _mic(centre - xyz[:, partner[0], :], box), axis=1
                        )
                        key = (f"{column}|{tag}|{c.partner_resname}"
                               f"{c.partner_canonical}:{c.partner_atom}")
                        per_frame[key] = d

            total_ns = _run_ns(Path(dcd).parent)
            if total_ns < args.min_ns:
                logging.warning("skipped %s %s rep%02d: only %.0f ns on disk",
                                mutation, state, rep, total_ns)
                continue
            mask = _window_mask(traj.n_frames, total_ns, args.window_ns)
            for f in np.flatnonzero(mask):
                f = int(f)
                rec = dict(mutation=mutation, state=state, replicate=rep, frame=f,
                           run_ns=round(total_ns, 1))
                for key, series in per_frame.items():
                    rec[key] = float(series[f])
                rows.append(rec)
            logging.info("%s %s rep%02d: %.0f ns on disk, %d of %d frames in the "
                         "last %.0f ns, %d surrogates", mutation, state, rep, total_ns,
                         int(mask.sum()), traj.n_frames, args.window_ns, len(per_frame))
        except Exception as exc:
            logging.error("FAILED %s %s rep%s: %s", mutation, state, rep, exc)

    if not rows:
        return 1
    df = pd.DataFrame(rows)
    suffix = "" if args.state == "holo" else f"_{args.state}"
    out = args.output_dir / f"tables/surrogate_per_frame{suffix}.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, index=False)
    logging.info("wrote %s (%d rows, %d cols)", out, len(df), df.shape[1])
    return 0


def _genotype_means(df: pd.DataFrame, cols: list[str]) -> pd.DataFrame:
    """Mean per replicate, then mean across replicates -- as Table 3 is built."""
    per_rep = df.groupby(["mutation", "replicate"])[cols].mean().reset_index()
    return per_rep.groupby("mutation")[cols].mean()


def stage_rank(args) -> int:
    from scipy.stats import spearmanr

    src = args.output_dir / "tables/surrogate_per_frame_both.csv"
    if not src.exists():
        src = args.output_dir / "tables/surrogate_per_frame.csv"
    sur = pd.read_csv(src)
    if "state" in sur.columns:
        sur = sur[sur.state == "holo"]      # Table 3 is defined in the holo state only
    mech = pd.read_csv(args.mechanisms)
    keys = ["mutation", "replicate", "frame"]
    coord_cols = sorted({co for _, co, _, _, _ in AXES})
    merged = sur.merge(mech[keys + coord_cols], on=keys, how="inner")
    logging.info("merged %d frames across %d genotypes",
                 len(merged), merged.mutation.nunique())

    sur_cols = [c for c in sur.columns if "|" in c]
    means = _genotype_means(merged, sur_cols + coord_cols)

    rng = np.random.default_rng(0)
    rows = []
    by_column: dict[str, list[str]] = {}
    for c in sur_cols:
        by_column.setdefault(c.split("|")[0], []).append(c)

    for column, coord, _, _, _ in AXES:
        cands = by_column.get(column, [])
        if not cands:
            continue
        target = means[coord]
        ok = target.notna()
        rhos = {}
        for c in cands:
            sub = ok & means[c].notna()
            if sub.sum() < 8:
                continue
            rho, pval = spearmanr(means.loc[sub, c], target[sub])
            # secondary: does it track within a trajectory? thermal coupling,
            # reported but never used to choose.
            within = []
            for (_, _), g in merged.groupby(["mutation", "replicate"]):
                if g[c].notna().sum() > 20 and g[coord].notna().sum() > 20:
                    wr, _ = spearmanr(g[c], g[coord], nan_policy="omit")
                    if np.isfinite(wr):
                        within.append(wr)
            rhos[c] = rho
            _, tag, partner = c.split("|")
            rows.append(dict(
                table3_column=column, coordinate=coord, anchor=tag, partner=partner,
                n_genotypes=int(sub.sum()),
                rho_across_genotypes=round(float(rho), 3),
                p_naive=round(float(pval), 4),
                rho_within_replicate_median=(round(float(np.median(within)), 3)
                                             if within else None),
                surrogate=c,
            ))

        # Selection-corrected p for the winner: how often does shuffling the
        # genotype labels produce a best-of-this-set as strong by chance?
        if rhos:
            obs = max(abs(v) for v in rhos.values())
            mat = means.loc[ok, list(rhos)].to_numpy()
            tgt = target[ok].to_numpy()
            null = np.empty(args.permutations)
            for i in range(args.permutations):
                perm = rng.permutation(tgt)
                null[i] = max(
                    abs(spearmanr(mat[:, j], perm).statistic)
                    for j in range(mat.shape[1])
                )
            p_sel = float((null >= obs).mean())
            best = max(rhos, key=lambda k: abs(rhos[k]))
            for r in rows:
                if r["surrogate"] == best:
                    r["is_best_for_column"] = True
                    r["p_selection_corrected"] = round(p_sel, 4)
            logging.info("%s: best |rho|=%.3f, selection-corrected p=%.4f (%s)",
                         column, obs, p_sel, best.split("|", 1)[1])

    out = pd.DataFrame(rows)
    if out.empty:
        logging.error("no surrogate/coordinate pair had enough genotypes in common; "
                      "is surrogate_per_frame.csv the full scan?")
        return 1
    if "is_best_for_column" in out.columns:
        out["is_best_for_column"] = out["is_best_for_column"].fillna(False)
    out = out.sort_values(["table3_column", "rho_across_genotypes"],
                          key=lambda s: s.abs() if s.name == "rho_across_genotypes" else s,
                          ascending=[True, False])
    dest = args.output_dir / "tables/surrogate_ranking.csv"
    out.to_csv(dest, index=False)
    logging.info("wrote %s (%d rows)", dest, len(out))

    for column, _, _, _, _ in AXES:
        sub = out[out.table3_column == column].head(4)
        if sub.empty:
            continue
        print(f"\n=== {column} ===")
        print(sub[["anchor", "partner", "n_genotypes", "rho_across_genotypes",
                   "p_naive", "rho_within_replicate_median",
                   "p_selection_corrected"]].to_string(index=False))
    return 0


def _group_stats(df: pd.DataFrame, col: str) -> dict:
    """Mean, SEM and replicate count for one (genotype, state) cell."""
    per_rep = df.groupby("replicate")[col].mean().dropna()
    n = len(per_rep)
    if n == 0:
        return dict(mean=np.nan, var=np.nan, n=0)
    sd = float(per_rep.std(ddof=1)) if n > 1 else np.nan
    return dict(mean=float(per_rep.mean()),
                var=(sd ** 2 / n) if n > 1 else np.nan, n=n)


def _welch(parts: list[dict]) -> tuple[float, float]:
    """Combined standard error and Welch-Satterthwaite df for a sum of means.

    Three replicates per cell means two degrees of freedom, so a 95% interval
    needs t = 4.30, not 1.96. Reporting +/- 2 SEM here would understate the
    band by more than a factor of two and make cells of the 2x2 look
    distinguishable when they are not.
    """
    varis = [p["var"] for p in parts]
    if any(not np.isfinite(v) for v in varis):
        return np.nan, np.nan
    total = float(sum(varis))
    denom = sum(v ** 2 / max(1, p["n"] - 1) for v, p in zip(varis, parts))
    df = (total ** 2 / denom) if denom > 0 else np.nan
    return float(np.sqrt(total)), float(df)


def stage_ddf(args) -> int:
    from scipy.stats import t as student_t

    src = args.output_dir / "tables/surrogate_per_frame_both.csv"
    df = pd.read_csv(src)
    chosen = pd.read_csv(args.surrogates) if args.surrogates else None
    if chosen is None:
        rank = pd.read_csv(args.output_dir / "tables/surrogate_ranking.csv")
        chosen = (rank[rank.get("is_best_for_column", False) == True]
                  if "is_best_for_column" in rank.columns else rank)
    features = list(dict.fromkeys(chosen.surrogate))
    logging.info("%d features, %d genotypes in scan", len(features), df.mutation.nunique())

    ddg = pd.read_csv(args.ddg).set_index("mutation")
    rows = []
    for feature in features:
        wt = {st: _group_stats(df[(df.mutation == "WT") & (df.state == st)], feature)
              for st in ("apo", "holo")}
        for mutation in sorted(set(df.mutation) - {"WT"}):
            mu = {st: _group_stats(df[(df.mutation == mutation) & (df.state == st)], feature)
                  for st in ("apo", "holo")}
            if any(c["n"] == 0 for c in list(mu.values()) + list(wt.values())):
                continue
            d_holo = mu["holo"]["mean"] - wt["holo"]["mean"]
            d_apo = mu["apo"]["mean"] - wt["apo"]["mean"]
            ddf = d_holo - d_apo
            se, dof = _welch([mu["holo"], wt["holo"], mu["apo"], wt["apo"]])
            crit = float(student_t.ppf(0.975, dof)) if np.isfinite(dof) and dof > 0 else np.nan
            ci = crit * se if np.isfinite(crit) else np.nan

            g = ddg.loc[mutation] if mutation in ddg.index else None
            ddg_mean = float(g["ddg_mean"]) if g is not None else np.nan
            ddg_sem = float(g["ddg_sem"]) if g is not None else np.nan
            ddg_n = int(g["n_replicates"]) if g is not None else 0
            ddg_ci = (float(student_t.ppf(0.975, ddg_n - 1)) * ddg_sem
                      if ddg_n > 1 and np.isfinite(ddg_sem) else np.nan)

            ddf_sig = bool(np.isfinite(ci) and abs(ddf) > ci)
            ddg_sig = bool(np.isfinite(ddg_ci) and abs(ddg_mean) > ddg_ci)
            if not np.isfinite(ci) or not np.isfinite(ddg_ci):
                cell = "undetermined"
            elif ddf_sig and not ddg_sig:
                cell = "3 differential, energy blind"
            elif not ddf_sig and not ddg_sig:
                cell = "1 cancels (or underpowered)"
            elif ddf_sig and ddg_sig:
                cell = "2 differential, energy sees it"
            else:
                cell = "4 feature misses it"

            rows.append(dict(
                feature=feature, mutation=mutation,
                fold_reduction=(float(g["fold_reduction"])
                                if g is not None and pd.notna(g["fold_reduction"]) else np.nan),
                n_apo=mu["apo"]["n"], n_holo=mu["holo"]["n"],
                delta_holo=round(d_holo, 3), delta_apo=round(d_apo, 3),
                ddf=round(ddf, 3), ddf_ci95=round(ci, 3) if np.isfinite(ci) else np.nan,
                ddf_significant=ddf_sig,
                ddg=round(ddg_mean, 3), ddg_ci95=round(ddg_ci, 3) if np.isfinite(ddg_ci) else np.nan,
                ddg_significant=ddg_sig,
                mismatch_genotype=mutation in MISMATCH, cell=cell,
            ))

    out = pd.DataFrame(rows)
    if out.empty:
        logging.error("no genotype had both states")
        return 1
    dest = args.output_dir / "tables/ddf_vs_ddg.csv"
    out.to_csv(dest, index=False)
    logging.info("wrote %s (%d rows)", dest, len(out))

    for feature, g in out.groupby("feature", sort=False):
        print(f"\n=== {feature.split('|',1)[1]}  (for {feature.split('|')[0]}) ===")
        cols = ["mutation", "fold_reduction", "delta_holo", "delta_apo",
                "ddf", "ddf_ci95", "ddg", "ddg_ci95", "mismatch_genotype"]
        print(g.sort_values("ddf")[cols].to_string(index=False))
        print("  cells:", dict(g.cell.value_counts()))
    return 0


def stage_plot(args) -> int:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    d = pd.read_csv(args.output_dir / "tables/ddf_vs_ddg.csv")
    feats = list(dict.fromkeys(d.feature))
    ncol = min(2, len(feats))
    nrow = int(np.ceil(len(feats) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(7.0 * ncol, 6.0 * nrow), squeeze=False)

    for ax, feature in zip(axes.ravel(), feats):
        g = d[d.feature == feature]
        # The bands are where a point cannot be told from zero, so the four
        # quadrants of the 2x2 are only meaningful outside them.
        xb = float(np.nanmedian(g.ddf_ci95))
        yb = float(np.nanmedian(g.ddg_ci95))
        ax.axvspan(-xb, xb, color="#b0b0b0", alpha=0.16, zorder=0)
        ax.axhspan(-yb, yb, color="#b0b0b0", alpha=0.16, zorder=0)
        ax.axvline(0, color="#555", lw=0.9, zorder=1)
        ax.axhline(0, color="#555", lw=0.9, zorder=1)

        for _, r in g.iterrows():
            mism = bool(r.mismatch_genotype)
            colour = "#c92a2a" if mism else "#2c7fb8"
            ax.errorbar(r.ddf, r.ddg, xerr=r.ddf_ci95, yerr=r.ddg_ci95,
                        fmt="o", ms=9 if mism else 6,
                        color=colour, ecolor=colour, elinewidth=1.0,
                        capsize=2.5, alpha=0.95 if mism else 0.65, zorder=4 if mism else 3)
            if mism or abs(r.ddf) > xb:
                ax.annotate(r.mutation, (r.ddf, r.ddg), xytext=(6, 5),
                            textcoords="offset points", fontsize=8,
                            color=colour, fontweight="bold" if mism else "normal")

        label = feature.split("|", 1)[1]
        ax.set_xlabel(f"$\\Delta\\Delta$ {label}   (holo $-$ apo, Å)")
        ax.set_ylabel(r"$\Delta\Delta G_{bind}$ (kcal/mol)")
        ax.set_title(f"{feature.split('|')[0]}  via  {label}", fontsize=10)
        ax.spines[["top", "right"]].set_visible(False)

    for ax in axes.ravel()[len(feats):]:
        ax.set_visible(False)

    from matplotlib.lines import Line2D
    fig.legend(handles=[
        Line2D([], [], marker="o", ls="none", color="#c92a2a",
               label="High resistance, $\\Delta\\Delta G_{bind}$ does not reproduce it"),
        Line2D([], [], marker="o", ls="none", color="#2c7fb8", label="Other genotypes"),
        matplotlib.patches.Patch(facecolor="#b0b0b0", alpha=0.3,
                                 label="Indistinguishable from zero (95%)"),
    ], loc="lower center", ncol=3, fontsize=9, frameon=False)
    fig.suptitle("Does the mutation reshape the pocket differently with and without "
                 "doravirine, and does the binding energy see it?", fontsize=12)
    fig.tight_layout(rect=(0, 0.05, 1, 0.97))
    dest = args.output_dir / "ddf_vs_ddg.png"
    fig.savefig(dest, dpi=300)
    fig.savefig(dest.with_suffix(".pdf"))
    logging.info("wrote %s", dest)
    return 0


def stage_axes(args) -> int:
    import mdtraj as md

    mf = pd.read_csv(args.manifest)
    wt = mf[(mf["mutation"] == "WT") & (mf["replicate"] == 1)]
    if wt.empty:
        logging.error("no WT replicate 1 in %s", args.manifest)
        return 1
    topo, dcd, _ = _replicate_inputs(wt.iloc[0], REPO_ROOT)
    traj = md.load(str(dcd), top=str(topo), stride=args.stride)
    traj.make_molecules_whole(inplace=True)
    rings = _ligand_rings(traj)

    out = []
    for column, coord, canonical, names, lig_key in AXES:
        anchor_idx = _anchor(traj.topology, canonical, names)
        lig_idx = rings.get(lig_key) or []
        if not anchor_idx or not lig_idx:
            logging.warning("%s: anchor or ligand group empty, skipped", column)
            continue
        cand = _candidates(traj, anchor_idx, lig_idx, canonical,
                           backbone_only=args.partner_atoms == "backbone")
        cand.insert(0, "table3_column", column)
        cand.insert(1, "coordinate", coord)
        out.append(cand)
        logging.info("%s: %d candidate partner atoms", column, len(cand))

    if not out:
        return 1
    df = pd.concat(out, ignore_index=True)
    args.candidates.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.candidates, index=False)
    logging.info("wrote %s (%d rows)", args.candidates, len(df))
    return 0


def main(argv: list[str] | None = None) -> int:
    warnings.filterwarnings("ignore")
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--stage", choices=("axes", "scan", "rank", "ddf", "plot"), default="axes")
    p.add_argument("--manifest", type=Path, default=Path("manifests/md_manifest.csv"))
    p.add_argument("--stride", type=int, default=1)
    p.add_argument("--mutations", nargs="*", default=None)
    p.add_argument("--mechanisms", type=Path,
                   default=Path("results/analysis/mechanisms/mechanism_coordinates.csv"))
    p.add_argument("--permutations", type=int, default=2000)
    p.add_argument("--window-ns", type=float, default=50.0,
                   help="Average over the final N ns of each run, so runs of "
                        "different length and stride are compared on equal footing.")
    p.add_argument("--min-ns", type=float, default=80.0,
                   help="Skip replicates with less than this much trajectory.")
    p.add_argument("--ddg", type=Path,
                   default=Path("results/analysis/binding_energy/tables/"
                                "mutation_ddg_summary.csv"))
    p.add_argument("--surrogates", type=Path, default=None,
                   help="CSV with a 'surrogate' column naming the features to test; "
                        "defaults to the winner per Table 3 column.")
    p.add_argument("--state", choices=("holo", "apo", "both"), default="holo",
                   help="Which simulations to scan. The surrogates are "
                        "protein-protein distances, so apo needs no ligand.")
    p.add_argument("--partner-atoms", choices=("backbone", "all"), default="backbone",
                   help="Restrict the far-wall partner atom to the backbone, so a "
                        "mutation at that position cannot move the surrogate by "
                        "side-chain identity alone.")
    p.add_argument("--output-dir", type=Path,
                   default=REPO_ROOT / "results/analysis/pocket_surrogates")
    args = p.parse_args(argv)
    args.candidates = args.output_dir / "tables/surrogate_candidates.csv"

    if args.stage == "axes":
        return stage_axes(args)
    if args.stage == "scan":
        return stage_scan(args)
    if args.stage == "rank":
        return stage_rank(args)
    if args.stage == "ddf":
        return stage_ddf(args)
    if args.stage == "plot":
        return stage_plot(args)
    logging.error("stage %s not implemented yet", args.stage)
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
