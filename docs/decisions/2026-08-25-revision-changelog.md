# Manuscript revision changelog

Source: `paper/submission/DorDRM-MD-08-14-26.docx`
Output: `paper/submission/DorDRM-MD-08-14-26-rev.docx`

- Rewrote Abstract (no regression claim; FEP incomplete; occupancy stats).
- Filled Intro FEP-methods paragraph; removed placeholder heading.
- Fixed Non-Nucleotide → non-nucleoside.
- Filled Table 1 category for K103N+M230L → Resistant.
- Demoted MM/GBSA heading; scrubbed FEP language from that block.
- Inserted FEP Results + NNIBP descriptor suite sections before contact analysis.
- Added replicate CI / Welch / descriptive-only language to occupancy screen.
- Downgraded Y188L Lys102/Pro225 occupancy claim to descriptive_only.
- Hedged V106A/G190E causal verbs; kept Ser105/Val179 as supported reporters.
- Marked V106I+F227C–227 occupancy as descriptive_only.
- Deleted linear-regression Results/Methods/Discussion blocks (13 paragraphs).
- Replaced residual regression Discussion with explicit no-model statement.
- Expanded Methods: FF justification, FEP protocol, occupancy stats, descriptor suite.
- Discussion lead reframed around incomplete FEP + local MD claims.

## Occupancy verdicts (for in-text numbers)

See `results/analysis/occupancy_stats/OCCUPANCY_STATS_NOTES.md`.

| callout | Δocc | verdict |
|---|---:|---|
| V106A / SER105 | +0.75 | supported |
| V106A combos / SER105 | +0.75–0.88 | supported |
| G190E / VAL179 | −0.84 | supported |
| V106A family / LYS104 | +0.39–0.67 | mixed (Welch only) |
| Y188L Lys102 / Pro225 | +0.19 / −0.13 | descriptive_only |
| V106I+F227C / 227 | −0.26 | descriptive_only |

## Still manual / figure swaps

- Drop Figure 8 (regression) from the Word gallery and renumber.
- Insert FEP protocol panels (`results/analysis/fep_pmx/protocol/V106A/`) and no-fit scatter (`panel_ddg_vs_experiment.png`).
- Insert occupancy CI plots (`results/analysis/occupancy_stats/plots/`).
- Insert modern_md_suite highlights (pose PCA, H-bonds, COM r–α, DCCM).
- Collapse Table 2 MM/GBSA to SI if desired; main text already demotes it.
- Add Gapsys/Aldeghi/Hauser bibliography entries if not already cited numerically.
- Re-run `python ops/maintenance/manuscript/apply_atanu_revision.py` only from the pristine `DorDRM-MD-08-14-26.docx` (script is not idempotent on `-rev.docx`).
