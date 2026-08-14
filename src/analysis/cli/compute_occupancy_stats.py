#!/usr/bin/env python3
"""Replicate CIs + Welch tests for highlighted DOR contact Δoccupancy.

Uses the existing residue-level replicate_contact table + global permutation null.
Downgrades text-highlighted shifts that fail both the global FWER null and the
per-contact Welch test (n=3 vs 3) to descriptive_only.

    python -m src.analysis.cli.compute_occupancy_stats
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

REPO = Path(__file__).resolve().parents[3]
DEFAULT_REP = Path(
    "results/analysis/triplet_story_analyses/"
    "contact_story_all_mutations_excluding_f227c/tables/replicate_contact.csv"
)
DEFAULT_PERM = Path(
    "results/analysis/triplet_story_analyses/"
    "contact_occupancy_permutation_null/tables/residue_shift_permutation_results.csv"
)
DEFAULT_OUT = Path("results/analysis/occupancy_stats")

# Plan §5 + draft story callouts. Keys are auth residue numbers.
HIGHLIGHTS: list[tuple[str, int, str, str]] = [
    # (mutation, auth_resid, residue_label, why)
    ("V106A", 105, "SER105", "V106A entrance-shift reporter"),
    ("V106A", 104, "LYS104", "V106A entrance-shift companion"),
    ("V106A+F227L", 105, "SER105", "V106A combo Ser105"),
    ("V106A+F227L", 104, "LYS104", "V106A combo Lys104"),
    ("V106A+L234I", 105, "SER105", "V106A combo Ser105"),
    ("V106A+L234I", 104, "LYS104", "V106A combo Lys104"),
    ("V106A+P225H", 105, "SER105", "V106A combo Ser105"),
    ("V106A+P225H", 104, "LYS104", "V106A combo Lys104"),
    ("V106I", 105, "SER105", "negative control vs V106A"),
    ("V106I+F227C", 227, "RES227", "227-region story"),
    ("V106I+F227C", 105, "SER105", "contrast vs V106A+F227L"),
    ("V106A+F227L", 227, "RES227", "contrast vs V106I+F227C"),
    ("Y188L", 102, "LYS102", "draft-highlighted Δoccupancy"),
    ("Y188L", 225, "PRO225", "draft-highlighted Δoccupancy"),
    ("Y188L", 188, "TYR188", "direct mutation site packing"),
    ("K103N+P225H", 225, "PRO225", "draft-highlighted Δoccupancy"),
    ("G190E", 179, "VAL179", "G190E hairpin reporter"),
    ("G190A", 179, "VAL179", "negative control vs G190E"),
    ("G190S", 179, "VAL179", "intermediate G190 comparison"),
]


def _ci95(vals: np.ndarray) -> tuple[float, float, float]:
    """mean, lo, hi for n=3 (Student-t). Falls back to mean±inf if n<2."""
    x = np.asarray(vals, dtype=float)
    x = x[np.isfinite(x)]
    n = len(x)
    mean = float(np.mean(x)) if n else float("nan")
    if n < 2:
        return mean, float("nan"), float("nan")
    sem = float(np.std(x, ddof=1) / np.sqrt(n))
    tcrit = float(stats.t.ppf(0.975, df=n - 1))
    return mean, mean - tcrit * sem, mean + tcrit * sem


def _welch(a: np.ndarray, b: np.ndarray) -> tuple[float, float]:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    if len(a) < 2 or len(b) < 2:
        return float("nan"), float("nan")
    res = stats.ttest_ind(a, b, equal_var=False)
    return float(res.statistic), float(res.pvalue)


def _rep_occupancy_auth(rep: pd.DataFrame, mutation: str, auth: int) -> np.ndarray:
    """Per-replicate occupancy for auth residue; missing contacts → 0 (heatmap convention)."""
    mut_reps = sorted(rep.loc[rep["mutation"] == mutation, "replicate"].dropna().unique().tolist())
    if not mut_reps:
        return np.array([], dtype=float)
    sub = rep[(rep["mutation"] == mutation) & (rep["auth_resid"] == auth)]
    by_rep = sub.groupby("replicate", as_index=True)["occupancy"].mean()
    return np.array([float(by_rep.get(r, 0.0)) for r in mut_reps], dtype=float)


def compute_highlight_stats(
    rep: pd.DataFrame,
    perm: pd.DataFrame,
    highlights: list[tuple[str, int, str, str]],
) -> pd.DataFrame:
    wt_cache: dict[int, np.ndarray] = {}
    rows = []
    for mut, auth, label, why in highlights:
        if auth not in wt_cache:
            wt_cache[auth] = _rep_occupancy_auth(rep, "WT", auth)
        wt = wt_cache[auth]
        mut_occ = _rep_occupancy_auth(rep, mut, auth)

        m_mean, m_lo, m_hi = _ci95(mut_occ)
        w_mean, w_lo, w_hi = _ci95(wt)
        if mut_occ.size and np.isfinite(w_mean):
            d_mean, d_lo, d_hi = _ci95(mut_occ - w_mean)
            delta_point = float(m_mean - w_mean)
        else:
            d_mean = d_lo = d_hi = delta_point = float("nan")
        t_stat, t_p = _welch(mut_occ, wt)

        pr = perm[(perm["mutation"] == mut) & (perm["auth_resid"] == auth)]
        if pr.empty:
            g_p = float("nan")
            g_shift = float("nan")
            g_pass = False
        else:
            g_p = float(pr.iloc[0]["global_fwer_p_for_cell"])
            g_shift = float(pr.iloc[0]["wt_referenced_occupancy_shift"])
            g_pass = bool(pr.iloc[0]["exceeds_global_threshold"])

        welch_pass = bool(np.isfinite(t_p) and t_p < 0.05)
        if g_pass and welch_pass:
            verdict = "supported"
        elif g_pass or welch_pass:
            verdict = "mixed"
        else:
            verdict = "descriptive_only"

        rows.append(
            {
                "mutation": mut,
                "auth_resid": auth,
                "residue": label,
                "rationale": why,
                "n_mut_reps": int(mut_occ.size),
                "n_wt_reps": int(wt.size),
                "occ_mut_mean": m_mean,
                "occ_mut_ci95_lo": m_lo,
                "occ_mut_ci95_hi": m_hi,
                "occ_wt_mean": w_mean,
                "occ_wt_ci95_lo": w_lo,
                "occ_wt_ci95_hi": w_hi,
                "delta_vs_wt_mean": delta_point,
                "delta_vs_wt_ci95_lo": d_lo if np.isfinite(d_lo) else float("nan"),
                "delta_vs_wt_ci95_hi": d_hi if np.isfinite(d_hi) else float("nan"),
                "welch_t": t_stat,
                "welch_p": t_p,
                "welch_pass_p05": welch_pass,
                "perm_delta": g_shift,
                "perm_global_fwer_p": g_p,
                "perm_exceeds_global_90pct": g_pass,
                "verdict": verdict,
            }
        )
    return pd.DataFrame(rows)


def plot_highlight_bars(df: pd.DataFrame, out: Path) -> None:
    import matplotlib.pyplot as plt

    order = df.copy()
    order["label"] = order["mutation"] + " / " + order["residue"]
    # sort: supported first, then mixed, then descriptive; within by |delta|
    rank = {"supported": 0, "mixed": 1, "descriptive_only": 2}
    order["rk"] = order["verdict"].map(rank)
    order["abs_delta"] = order["delta_vs_wt_mean"].abs()
    order = order.sort_values(["rk", "abs_delta"], ascending=[True, False])

    colors = {
        "supported": "#c44e52",
        "mixed": "#dd8452",
        "descriptive_only": "#8da0cb",
    }
    fig, ax = plt.subplots(figsize=(10, max(5.5, 0.38 * len(order) + 1.2)))
    y = np.arange(len(order))
    x = order["delta_vs_wt_mean"].to_numpy()
    # error from CI half-width when available
    lo = order["delta_vs_wt_ci95_lo"].to_numpy()
    hi = order["delta_vs_wt_ci95_hi"].to_numpy()
    xerr = np.vstack([x - lo, hi - x])
    xerr = np.where(np.isfinite(xerr), xerr, 0.0)
    c = [colors[v] for v in order["verdict"]]
    ax.barh(y, x, xerr=xerr, color=c, ecolor="#333333", capsize=2, height=0.72, edgecolor="#333333", linewidth=0.4)
    ax.axvline(0, color="#333333", lw=1.0)
    ax.set_yticks(y, order["label"].tolist(), fontsize=9)
    ax.set_xlabel("Δ occupancy (mutant − WT mean) ± 95% CI (n = 3)")
    ax.set_title(
        "Highlighted DOR contact Δoccupancy with replicate CIs\n"
        "red = global-null + Welch; orange = one of two; blue = descriptive only",
        fontweight="bold",
        fontsize=11,
    )
    ax.grid(axis="x", alpha=0.25, linestyle=":")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(out, dpi=200)
    plt.close(fig)


def plot_story_ci_panels(rep: pd.DataFrame, df: pd.DataFrame, out: Path) -> None:
    """Story-focused Δoccupancy with CI whiskers (Fig-4 style update)."""
    import matplotlib.pyplot as plt

    stories = [
        ("Y188L", ["LYS102", "PRO225"], "Y188L"),
        ("V106A", ["SER105", "LYS104"], "V106A"),
        ("V106A+F227L", ["SER105", "LYS104"], "V106A+F227L"),
        ("V106I+F227C", ["RES227", "SER105"], "V106I+F227C"),
        ("G190E", ["VAL179"], "G190E"),
        ("K103N+P225H", ["PRO225"], "K103N+P225H"),
    ]
    fig, axes = plt.subplots(2, 3, figsize=(11, 7), squeeze=False)
    for ax, (mut, residues, title) in zip(axes.ravel(), stories):
        sub = df[(df["mutation"] == mut) & (df["residue"].isin(residues))].copy()
        if sub.empty:
            ax.set_title(title)
            ax.text(0.5, 0.5, "no data", ha="center", va="center", transform=ax.transAxes)
            continue
        sub = sub.set_index("residue").reindex(residues).reset_index()
        x = np.arange(len(sub))
        y = sub["delta_vs_wt_mean"].to_numpy()
        lo = sub["delta_vs_wt_ci95_lo"].to_numpy()
        hi = sub["delta_vs_wt_ci95_hi"].to_numpy()
        yerr = np.vstack([y - lo, hi - y])
        yerr = np.where(np.isfinite(yerr), yerr, 0.0)
        cols = [
            "#c44e52" if v == "supported" else "#dd8452" if v == "mixed" else "#8da0cb"
            for v in sub["verdict"]
        ]
        ax.bar(x, y, yerr=yerr, color=cols, ecolor="#333", capsize=3, width=0.65)
        ax.axhline(0, color="#333", lw=0.8)
        ax.set_xticks(x, residues, fontsize=9)
        ax.set_ylabel("Δocc ± 95% CI")
        ax.set_title(title, fontweight="bold")
        ax.grid(axis="y", alpha=0.25, linestyle=":")
        # annotate verdict letter
        for i, (_, row) in enumerate(sub.iterrows()):
            v = row["verdict"]
            if not isinstance(v, str) or v not in ("supported", "mixed", "descriptive_only"):
                continue
            tag = {"supported": "S", "mixed": "M", "descriptive_only": "D"}[v]
            yi = float(y[i]) if np.isfinite(y[i]) else 0.0
            ax.text(i, yi + (0.02 if yi >= 0 else -0.02), tag, ha="center", va="bottom" if yi >= 0 else "top", fontsize=8)
    fig.suptitle(
        "Story Δoccupancy with replicate 95% CIs  ·  S=supported  M=mixed  D=descriptive only",
        fontweight="bold",
    )
    fig.tight_layout()
    fig.savefig(out, dpi=200)
    plt.close(fig)


def write_notes(df: pd.DataFrame, out: Path) -> None:
    lines = [
        "# Occupancy statistics (Atanu #5)",
        "",
        "Replicate mean ± 95% CI (Student-t, n=3) and Welch t-test vs WT for draft-highlighted",
        "residue contacts. Global FWER p from the existing trajectory-label permutation null.",
        "",
        "## Verdict key",
        "",
        "- **supported** — exceeds global 90th-percentile null **and** Welch p < 0.05",
        "- **mixed** — passes exactly one of the two",
        "- **descriptive_only** — fails both; do not treat as significant in main text",
        "",
        "## Callouts for draft",
        "",
    ]
    for verdict in ("supported", "mixed", "descriptive_only"):
        sub = df[df["verdict"] == verdict].sort_values("delta_vs_wt_mean", key=lambda s: s.abs(), ascending=False)
        lines.append(f"### {verdict}")
        lines.append("")
        if sub.empty:
            lines.append("(none)")
            lines.append("")
            continue
        for _, r in sub.iterrows():
            lines.append(
                f"- **{r['mutation']} / {r['residue']}**: Δ = {r['delta_vs_wt_mean']:+.3f} "
                f"[{r['delta_vs_wt_ci95_lo']:+.3f}, {r['delta_vs_wt_ci95_hi']:+.3f}]; "
                f"Welch p = {r['welch_p']:.3g}; global FWER p = {r['perm_global_fwer_p']:.3g}"
            )
        lines.append("")
    lines.extend(
        [
            "## Draft language guidance",
            "",
            "- Keep Ser105 (V106A family) and Val179 (G190E) as load-bearing reporters if supported/mixed.",
            "- Downgrade Y188L Lys102 / Pro225 and V106I+F227C Phe227 if descriptive_only — "
            "phrase as ‘observed shift, not significant under replicate tests’.",
            "- n = 3 CIs are wide; never write ‘p confirms mechanism’.",
            "",
        ]
    )
    out.write_text("\n".join(lines))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--replicate-contact-csv", type=Path, default=DEFAULT_REP)
    ap.add_argument("--permutation-csv", type=Path, default=DEFAULT_PERM)
    ap.add_argument("--output-dir", type=Path, default=DEFAULT_OUT)
    args = ap.parse_args()
    out = args.output_dir
    (out / "tables").mkdir(parents=True, exist_ok=True)
    (out / "plots").mkdir(parents=True, exist_ok=True)

    rep = pd.read_csv(args.replicate_contact_csv)
    rep["mutation"] = rep["mutation"].astype(str)
    rep["auth_resid"] = pd.to_numeric(rep["auth_resid"], errors="coerce").astype("Int64")
    rep["replicate"] = pd.to_numeric(rep["replicate"], errors="coerce").astype("Int64")
    rep["occupancy"] = pd.to_numeric(rep["occupancy"], errors="coerce")

    perm = pd.read_csv(args.permutation_csv)
    perm["auth_resid"] = pd.to_numeric(perm["auth_resid"], errors="coerce").astype(int)
    df = compute_highlight_stats(rep, perm, HIGHLIGHTS)
    df.to_csv(out / "tables" / "highlighted_occupancy_stats.csv", index=False)

    # also dump full per-rep table for highlighted pairs
    rows = []
    for mut, auth, label, why in HIGHLIGHTS:
        for genotype in ("WT", mut):
            vals = _rep_occupancy_auth(rep, genotype, auth)
            for i, v in enumerate(vals, start=1):
                rows.append(
                    {
                        "mutation": genotype,
                        "auth_resid": auth,
                        "residue": label,
                        "replicate": i,
                        "occupancy": float(v),
                        "pair_mut": mut,
                        "rationale": why,
                    }
                )
    pd.DataFrame(rows).to_csv(out / "tables" / "highlighted_occupancy_per_rep.csv", index=False)

    plot_highlight_bars(df, out / "plots" / "highlighted_delta_occupancy_ci.png")
    plot_story_ci_panels(rep, df, out / "plots" / "story_delta_occupancy_ci.png")
    write_notes(df, out / "OCCUPANCY_STATS_NOTES.md")

    (out / "config" / "run_config.json").parent.mkdir(parents=True, exist_ok=True)
    (out / "config" / "run_config.json").write_text(
        json.dumps(
            {
                "replicate_contact_csv": str(args.replicate_contact_csv),
                "permutation_csv": str(args.permutation_csv),
                "highlights": [
                    {"mutation": m, "auth_resid": a, "residue": lab, "rationale": w}
                    for m, a, lab, w in HIGHLIGHTS
                ],
                "n_highlights": len(HIGHLIGHTS),
            },
            indent=2,
        )
    )
    print(df.groupby("verdict").size().to_string())
    print(f"Wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
