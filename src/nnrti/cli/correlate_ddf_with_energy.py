#!/usr/bin/env python3
"""Does binding-coupled remodeling of the NNIBP track the calculated energetics?

Reviewer 1 (comment 3) asks whether a mutation reshapes RT in the apo state as
well as the holo state, and whether such remodeling is therefore invisible to
ddG_bind. The thermodynamic cycle introduces each mutation twice, once in each
state, so whatever a mutation does equally to both ends cancels exactly. The
part that does not cancel is the double difference

    ddF = [F(mut, holo) - F(WT, holo)] - [F(mut, apo) - F(WT, apo)]

for a structural feature F measurable in both states. Seven of the nine Table 3
columns are protein-ligand contacts and have no apo counterpart, so each is
replaced by the intra-protein surrogate validated in Supplementary Figure S1;
V(NNIBP) is a protein-only quantity and is used directly.

ddF is then tested against both calculated energies: ddE_Total, the MM/GBSA
end-point estimate, and ddG_bind, the non-equilibrium alchemical free energy.
Four outcomes are possible per genotype, and they say different things:

    ddF ~ 0, dd(energy) ~ 0   the change cancels in the cycle -- the reviewer's
                              hypothesis, and the energy is blind by construction
    ddF != 0, dd(energy) != 0 binding-coupled remodeling the energy registers
    ddF != 0, dd(energy) ~ 0  remodeling the energy does not see, which refutes
                              cancellation and points somewhere else: either the
                              energy model or a mechanism that is not affinity
    ddF ~ 0, dd(energy) != 0  this feature set does not contain the coordinate

A correlation across genotypes answers the same question panel-wide rather than
one genotype at a time. No correlation together with energies that are
themselves near zero is the informative case: it would mean resistance in this
panel is not carried by the pocket geometry these features describe.

Side-chain identity cancels in ddF. Replacing a side chain shifts the feature in
the apo and the holo state alike, so the double difference removes it; that is
why an anchor on a mutated residue is admissible here.

F227C is an alchemical intermediate rather than a clinical genotype and carries
no fold-resistance value, so it is excluded, as it is from Table 3.

Usage
-----
    PYTHONPATH=src python -m nnrti.cli.correlate_ddf_with_energy
"""
from __future__ import annotations

import argparse
import logging
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import t as student_t

from nnrti.cli.find_pocket_surrogates import MISMATCH, _group_stats, _welch
from nnrti.cli.plot_surrogate_validation import (
    AXIS_LABEL_FONTSIZE,
    COLORS,
    LEGEND_FONTSIZE,
    MARKERS,
    POINT_LABEL_FONTSIZE,
    STAT_FONTSIZE,
    TICK_LABEL_FONTSIZE,
    WT_COLOR,
    category_for,
    surrogate_label,
)
from nnrti.fep.combine_neq import _repel_labels

REPO_ROOT = Path(__file__).resolve().parents[3]

VOLUME = "V(NNIBP)"
EXCLUDE = {"F227C"}

#: (column in the merged table, axis label, short name for the stats table)
ENERGIES = [
    ("dde", r"$\Delta\Delta E_\mathrm{Total}$  (kcal mol$^{-1}$)", "ddE_Total"),
    ("ddg", r"$\Delta\Delta G_\mathrm{bind}$  (kcal mol$^{-1}$)", "ddG_bind"),
]


def _limits(ax, xs: np.ndarray, ys: np.ndarray, pad: float = 0.18) -> None:
    """Axis limits from the point cloud, always including the origin."""
    for lo, hi, setter in ((min(xs.min(), 0.0), max(xs.max(), 0.0), ax.set_xlim),
                           (min(ys.min(), 0.0), max(ys.max(), 0.0), ax.set_ylim)):
        span = hi - lo
        margin = pad * span if span > 0 else 1.0
        setter(lo - margin, hi + margin)


def _feature_label(row) -> str:
    if row.surrogate == "DIRECT":
        return r"$\Delta\Delta V$(NNIBP)  (Å$^3$)"
    return rf"$\Delta\Delta${surrogate_label(row.table3_column, row.partner)}  (Å)"


def build_ddf(frames: pd.DataFrame, features: list[str]) -> pd.DataFrame:
    """ddF with a Welch-Satterthwaite 95% interval, one row per genotype."""
    rows = []
    for feature in features:
        wt = {st: _group_stats(frames[(frames.mutation == "WT") & (frames.state == st)], feature)
              for st in ("apo", "holo")}
        for mutation in sorted(set(frames.mutation) - {"WT"} - EXCLUDE):
            mu = {st: _group_stats(
                      frames[(frames.mutation == mutation) & (frames.state == st)], feature)
                  for st in ("apo", "holo")}
            if any(c["n"] == 0 for c in list(mu.values()) + list(wt.values())):
                continue
            d_holo = mu["holo"]["mean"] - wt["holo"]["mean"]
            d_apo = mu["apo"]["mean"] - wt["apo"]["mean"]
            se, dof = _welch([mu["holo"], wt["holo"], mu["apo"], wt["apo"]])
            crit = float(student_t.ppf(0.975, dof)) if np.isfinite(dof) and dof > 0 else np.nan
            rows.append(dict(
                feature=feature, mutation=mutation,
                n_apo=mu["apo"]["n"], n_holo=mu["holo"]["n"],
                delta_holo=d_holo, delta_apo=d_apo, ddf=d_holo - d_apo,
                ddf_ci95=crit * se if np.isfinite(crit) else np.nan,
                mismatch_genotype=mutation in MISMATCH,
            ))
    return pd.DataFrame(rows)


def main(argv: list[str] | None = None) -> int:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    root = REPO_ROOT / "results/analysis/pocket_surrogates"
    ap.add_argument("--surrogate-frames", type=Path,
                    default=root / "tables/surrogate_per_frame_both.csv")
    ap.add_argument("--volume-frames", type=Path,
                    default=root / "tables/nnibp_volume_per_frame.csv")
    ap.add_argument("--selected", type=Path, default=root / "tables/surrogates_selected.csv")
    ap.add_argument("--fep", type=Path,
                    default=REPO_ROOT / "results/analysis/fep_pmx/panel_ddg.csv",
                    help="Non-equilibrium alchemical ddG_bind, as published in Table 2.")
    ap.add_argument("--mmgbsa", type=Path,
                    default=REPO_ROOT / "results/analysis/binding_energy/tables/"
                                        "mutation_ddg_summary.csv",
                    help="MM/GBSA ddE_Total, as published in Table 2.")
    ap.add_argument("--output", type=Path, default=root / "ddf_vs_energy.png")
    ap.add_argument("--dpi", type=int, default=300)
    args = ap.parse_args(argv)

    sel = pd.read_csv(args.selected)
    sur = pd.read_csv(args.surrogate_frames)
    ddf = build_ddf(sur, [s for s in sel.surrogate if s != "DIRECT"])

    if args.volume_frames.exists():
        vol = pd.read_csv(args.volume_frames).rename(columns={"pocket_volume": VOLUME})
        ddf = pd.concat([ddf, build_ddf(vol, [VOLUME])], ignore_index=True)
    else:
        logging.warning("%s missing; V(NNIBP) left out", args.volume_frames)

    # ddE_Total: the MM/GBSA sum, Table 2 column 3. Its CI uses the replicate
    # SEM over three replicates, so t = 4.30 as everywhere else here.
    mm = pd.read_csv(args.mmgbsa).set_index("mutation")
    fep = pd.read_csv(args.fep).set_index("genotype")

    def _ci(sem, n):
        return (float(student_t.ppf(0.975, n - 1)) * sem
                if n > 1 and np.isfinite(sem) else np.nan)

    ddf["dde"] = [mm.ddg_mean.get(m, np.nan) for m in ddf.mutation]
    ddf["dde_ci95"] = [_ci(mm.ddg_sem.get(m, np.nan), int(mm.n_replicates.get(m, 0)))
                       for m in ddf.mutation]
    ddf["ddg"] = [fep.ddg_bind_kcal.get(m, np.nan) for m in ddf.mutation]
    ddf["ddg_ci95"] = [_ci(fep.sem_kcal.get(m, np.nan), int(fep.n_reps.get(m, 0)))
                       for m in ddf.mutation]
    ddf["fold_reduction"] = [fep.dor_fold_reduction.get(m, np.nan) for m in ddf.mutation]

    dest = args.output.parent / "tables/ddf_vs_energy.csv"
    dest.parent.mkdir(parents=True, exist_ok=True)
    ddf.round(3).to_csv(dest, index=False)
    logging.info("wrote %s (%d rows)", dest, len(ddf))

    order = [s for s in sel.surrogate if s != "DIRECT"]
    if VOLUME in set(ddf.feature):
        order.append(VOLUME)
    label_of = {r.surrogate: _feature_label(r) for _, r in sel.iterrows()}
    label_of[VOLUME] = _feature_label(sel[sel.surrogate == "DIRECT"].iloc[0])

    nrow, ncol = len(order), len(ENERGIES)
    fig, axes = plt.subplots(nrow, ncol, figsize=(6.6 * ncol, 5.0 * nrow), squeeze=False)
    stats_rows = []

    for i, feature in enumerate(order):
        g = ddf[ddf.feature == feature]
        # Whether ddF ~ 0 means cancellation or only a wide interval is settled
        # by the apo/holo diagonal: regressing the holo shift on the apo shift,
        # a slope near 1 says the mutation moves the feature by the same amount
        # in both states, which is what makes the cycle blind to it.
        diag = stats.linregress(g.delta_apo.to_numpy(float), g.delta_holo.to_numpy(float))
        for j, (ecol, elabel, ename) in enumerate(ENERGIES):
            ax = axes[i][j]
            sub = g.dropna(subset=["ddf", ecol])
            xs = sub.ddf.to_numpy(float)
            ys = sub[ecol].to_numpy(float)
            names = list(sub.mutation)

            lr = stats.linregress(xs, ys)
            loo = []
            for k in range(len(xs)):
                keep = np.ones(len(xs), bool); keep[k] = False
                loo.append((stats.linregress(xs[keep], ys[keep]).rvalue ** 2, names[k]))
            loo_min, loo_min_g = min(loo)
            stats_rows.append(dict(
                feature=feature, energy=ename, n=len(xs),
                slope=round(lr.slope, 4), r2=round(lr.rvalue ** 2, 3), p_value=lr.pvalue,
                r2_min_leave_one_out=round(loo_min, 3), most_influential=loo_min_g,
                n_ddf_resolved=int((sub.ddf.abs() > sub.ddf_ci95).sum()),
                n_energy_resolved=int((sub[ecol].abs() > sub[f"{ecol}_ci95"]).sum()),
                median_abs_ddf=round(float(sub.ddf.abs().median()), 3),
                median_ddf_ci95=round(float(sub.ddf_ci95.median()), 3),
                apo_holo_slope=round(diag.slope, 3),
                apo_holo_r2=round(diag.rvalue ** 2, 3),
                apo_holo_p=diag.pvalue,
            ))

            grid = np.linspace(xs.min(), xs.max(), 50)
            ax.plot(grid, lr.intercept + lr.slope * grid, "-", color="0.35", lw=1.4, zorder=2)
            ax.axhline(0.0, color="0.85", lw=0.9, zorder=1)
            ax.axvline(0.0, color="0.85", lw=0.9, zorder=1)

            seen = set()
            for _, r in sub.iterrows():
                cat = category_for(r.mutation)
                ax.errorbar(r.ddf, r[ecol], xerr=r.ddf_ci95, yerr=r[f"{ecol}_ci95"],
                            fmt=MARKERS.get(cat, "o"),
                            color=WT_COLOR if cat == "Wild type" else COLORS[cat],
                            ecolor="0.78", elinewidth=0.9, capsize=2, ms=7, zorder=3,
                            clip_on=True,
                            label=cat if cat not in seen else None)
                seen.add(cat)

            # Scale the axes to the estimates, not to the intervals. A cell
            # built from four means of three replicates has about two degrees
            # of freedom, so the widest 95% bands here run to +/-15 A; letting
            # them set the limits would compress every point into the middle
            # and hide the very scatter the panel is about. The bands are still
            # drawn, and run off the axis where they are that wide.
            _limits(ax, xs, ys)
            _repel_labels(ax, fig, list(xs), list(ys), names,
                          fontsize=POINT_LABEL_FONTSIZE, n_iter=600)
            ax.set_xlabel(label_of[feature], fontsize=AXIS_LABEL_FONTSIZE)
            ax.set_ylabel(elabel, fontsize=AXIS_LABEL_FONTSIZE)
            ax.tick_params(axis="both", labelsize=TICK_LABEL_FONTSIZE)
            pstr = "$P$ < 0.001" if lr.pvalue < 1e-3 else f"$P$ = {lr.pvalue:.3f}"
            ax.annotate(f"$R^2$ = {lr.rvalue**2:.2f},  {pstr},  $n$ = {len(xs)}",
                        xy=(0.5, 1.015), xycoords="axes fraction", va="bottom",
                        ha="center", fontsize=STAT_FONTSIZE)
            for side in ("top", "right"):
                ax.spines[side].set_visible(False)

    handles, labels = [], []
    for a in axes.ravel():
        for h, l in zip(*a.get_legend_handles_labels()):
            if l not in labels:
                handles.append(h); labels.append(l)
    keep = [l for l in ("Susceptible", "Uncertain", "Resistant") if l in labels]
    fig.legend([handles[labels.index(l)] for l in keep], keep, loc="lower center",
               ncol=len(keep), fontsize=LEGEND_FONTSIZE, frameon=False,
               bbox_to_anchor=(0.5, 0.0))

    fig.tight_layout(rect=(0, 0.03, 1, 0.99), h_pad=2.4, w_pad=2.0)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=args.dpi)
    fig.savefig(args.output.with_suffix(".pdf"))
    table = pd.DataFrame(stats_rows)
    # One test per feature per energy, so a nominal P has to be paid for across
    # the whole grid. Holm rather than Bonferroni: it is uniformly the stronger
    # of the two and makes no independence assumption the tests do not meet.
    rank = table.p_value.rank(method="first").astype(int)
    adj = (len(table) - rank + 1) * table.p_value
    table["p_holm"] = np.minimum(
        1.0, adj.sort_index(key=lambda i: rank[i]).cummax().reindex(table.index))
    sdest = args.output.parent / "tables/ddf_vs_energy_stats.csv"
    table.to_csv(sdest, index=False)
    logging.info("wrote %s and %s", args.output, sdest)
    print(table.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
