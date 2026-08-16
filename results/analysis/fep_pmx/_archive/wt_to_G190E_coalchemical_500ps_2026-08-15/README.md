# Archived: wt_to_G190E, co-alchemical ion protocol, 500 ps switches

Archived 2026-08-15, before the G190E re-run under the K103N-matched protocol.
**Nothing here is a valid panel result. Do not feed it to `combine_neq`.**

## What this is

The G190E leg as run under the **abandoned co-alchemical ion** protocol
(`USE_COALCHEMICAL_ION = True`), with 500 ps switches:

| | value |
|---|---|
| ΔΔG_bind | +2.50 kcal/mol |
| SEM | 1.41 kcal/mol |
| switch_ps | 500.0 |
| per-rep apo BAR error | up to **35.49** kcal/mol (rep 1) |

That per-rep error is the signature of the non-convergence documented in
`scripts/fep_pmx/OPERATIONS.md` §7: decoupling a whole Cl- dissipates ~20-26
kcal/mol, driving forward/reverse work distributions ~20 kcal apart with
near-zero Crooks overlap, so BAR cannot converge.

## Why it was archived rather than deleted

`combine_neq` is now run **without `--force`** (the raw `dgdl.xvg` for the other
genotypes was destroyed by a `git clean -fd` on 2026-08-13, so `--force` would
fail on missing inputs). Without `--force` it reuses any cached `analysis.json`
it finds. Left in place, these files would have put a stale, non-convergent
G190E onto `panel_ddg.csv` and the panel scatter with no warning.

Kept because it is the empirical record of *why* the co-alchemical ion was
dropped, and it pairs with the K103N co-alchemical data as the before/after for
that decision.

## What was NOT archived (and the trap it set)

Only the per-unit `neq/analysis/` directories and `targets/G190E/summary.json`
were moved. The rest of `legs/wt_to_G190E/*/rep_*/neq/` — including
`neq_prepare.json` (`switch_ps: 500.0`) and the rendered 500 ps
`mdp/nonequil_{fwd,rev}.mdp` — was left in place on the assumption that a rebuild
would overwrite it. **It does not.** `prepare_neq.py` skips any unit that already
has a `neq_manifest.csv`, so the first prep run wrote only the panel manifest and
kept the 500 ps config. Re-running with `--force` fixes it.

If you archive another leg's results, move or delete the whole `neq/` subtree,
not just `analysis/`.

## What replaces it

G190E re-run in a **raw non-neutral box + Rocklin/Hunenberger analytical
net-charge correction** (`scripts/fep_pmx/charge_correction.py`), at **100 ps**
switches / 100 snapshots / 5 ns equilibration / 3 reps — identical to the
protocol K103N and the three K103N-compound legs ran. `wt_to_G190E` was removed
from `_BASE_LONG_SWITCH_LEGS` in `config.py` for that reason.
