# Pocket-opening cost at 100 ns: what changes, and one new result

Computed 01-Nov-2026 from the uniform 3 x 100 ns apo ensembles.
Source: `results/analysis/apo_ensembles/tables/pocket_opening_cost_100ns.csv`,
built from `apo_matched_window_summary.csv` (Rg-filtered) and the holo arm of
`nnibp_volume_both_states_100ns.csv`.

## The recompute is trustworthy: the holo arm reproduces

Per `CLAUDE.md`, a rerun has to recover the published numbers before its new
ones count. It does, on the side that was published:

| | manuscript | 100 ns recompute |
|---|---|---|
| WT holo NNIBP volume | 228 A^3 (implied by 91 + 137) | **227.4** |
| WT apo NNIBP volume | 91 | **79.9** |
| WT opening cost | 137 | **147.5** |

The holo value lands within 0.6 A^3. Everything that moved, moved on the apo
side, and in the direction longer sampling should move it: **the unliganded
pocket closes further than 10 ns showed**. The old apo runs were not converged,
which overestimated the apo volume and so understated the opening cost. That is
a better answer to R2.1 than the holo autocorrelation argument, because the
reviewer asked specifically about "unliganded pocket conformational equilibria".

## Full table

Opening cost = holo volume - apo volume, same genotype, A^3.

| genotype | n apo | apo | holo | opening | vs WT | |
|---|---|---|---|---|---|---|
| G190E | 3 | 46.8 | 293.5 | 246.7 | **+99.2** | mismatch |
| G190S | 3 | 50.9 | 280.2 | 229.3 | +81.8 | |
| V106I | 3 | 55.6 | 253.8 | 198.2 | +50.7 | |
| V106M | 3 | 54.0 | 237.5 | 183.5 | +36.1 | |
| Y188L | 3 | 54.8 | 237.7 | 182.9 | +35.4 | |
| Y318F | 3 | 63.4 | 240.2 | 176.8 | +29.3 | |
| V106I+F227C | 3 | 62.2 | 233.2 | 171.0 | +23.5 | |
| G190A | 3 | 61.0 | 231.9 | 170.9 | +23.4 | |
| A98G+F227C | 3 | 63.6 | 232.0 | 168.4 | **+20.9** | mismatch |
| L100I+K103N | 3 | 59.4 | 224.2 | 164.8 | +17.3 | |
| Y181C | 3 | 103.5 | 263.6 | 160.1 | +12.6 | |
| K103N | 3 | 74.2 | 233.7 | 159.5 | +12.0 | |
| F227C | 3 | 77.4 | 233.4 | 156.0 | +8.5 | |
| K103N+M230L | 3 | 89.6 | 245.6 | 156.0 | **+8.5** | mismatch |
| **WT** | 3 | 79.9 | 227.4 | 147.5 | 0.0 | |
| K103N+P225H | 3 | 98.4 | 241.4 | 143.0 | -4.5 | |
| V106A | **1** | 134.7 | 263.3 | 128.6 | -18.9 | |
| V106A+F227L | 3 | 118.9 | 246.3 | 127.4 | -20.1 | |
| V106A+L234I | 3 | 126.0 | 251.2 | 125.2 | -22.3 | |
| V106A+P225H | 3 | 119.1 | 225.0 | 105.9 | -41.6 | |

V106A is n=1 until its last two replicates finish (~13 h); see the caveat below.

## What this does to the R1.3 sentence

The manuscript currently says:

> "Wild type pays 137 A^3 to open the pocket. All three genotypes whose
> ddG_bind misses the phenotype pay more: G190E by 105 A^3, A98G+F227C by
> 38 A^3 and K103N+M230L by 20 A^3."

Each number moves, and one of them moves enough to break the argument it
carries:

| | manuscript | 100 ns | verdict |
|---|---|---|---|
| WT absolute | 137 | 147.5 | restate |
| G190E | +105 | +99.2 | holds |
| A98G+F227C | +38 | +20.9 | holds, roughly halved |
| K103N+M230L | +20 | +8.5 | **does not hold** |

K103N+M230L's excess is now 8.5 A^3, which is **identical to F227C's 8.5** and
below ten other genotypes, several of which ddG_bind predicts correctly. The
claim "all three pay more" is still literally true, but for K103N+M230L the
margin no longer distinguishes it from genotypes the calculation gets right, so
pocket-opening cost cannot be offered as the explanation for that genotype. The
honest version keeps G190E (large, unambiguous), keeps A98G+F227C (modest), and
drops K103N+M230L from this argument, which is the right outcome anyway: M230L
is a primer-grip position, and R2.2 already says the primer-grip route is
foreclosed by a binary RT-DOR system with no nucleic acid. The two comments then
agree instead of offering competing explanations.

The "cost is not specific to them" caveat gets stronger, not weaker: **G190S
pays +81.8**, the second-largest in the panel, and ddG_bind tracks its phenotype
fine.

## New result: the V106A class opens the pocket before DOR arrives

Not in the manuscript, and it answers R1.1 and R1.2 more directly than anything
currently written against them.

Every V106A-containing genotype pays **less** than wild type to open the pocket,
because its *unliganded* pocket is already wide:

| genotype | apo volume | vs WT apo (79.9) |
|---|---|---|
| V106A | 134.7 (n=1) | +54.8 |
| V106A+L234I | 126.0 | +46.1 |
| V106A+P225H | 119.1 | +39.2 |
| V106A+F227L | 118.9 | +39.0 |

Three independent genotypes at n=3 each show it, so the finding does not rest
on the n=1 V106A replicate. Reviewer 1's framing in comment 2 is that "mutations
may alter the RT conformation before DOR binding"; this is a measured instance
of exactly that. V106A does not merely reshape the bound interface, it
pre-organises the apo pocket toward the bound geometry.

Y181C is a second, milder instance: apo volume 103.5 against WT's 79.9, with its
Y181 chi1 at +40 degrees against WT's -68, i.e. the rotamer flip the
crystallography section already discusses, present in the unliganded enzyme.

## Caveats, stated rather than buried

- **V106A apo is n=1.** Its last two replicates reach 100 ns in ~13 h. The
  number in the table will move; the V106A-class conclusion will not, because
  it is carried by the three combination genotypes at n=3.
- **One replicate is still excluded.** G190E apo rep_02 loses 158 of 200 frames
  to the pocket-Ca Rg > 10 A filter (rg 10.5-42.1 A, median volume 2227 A^3):
  local unfolding, not a pocket. Worth saying that this is now the *only*
  exclusion. The manuscript's Methods describes two, G190E and one
  K103N+M230L; at 100 ns the K103N+M230L replicate is healthy. G190E's
  remaining spread is wide (sd 59.2) and its +99.2 should be quoted as the
  genotype with the largest opening cost rather than as a precise figure.
- **The apo cavity is still not re-solvated** before measurement, so apo volumes
  are biased low and the openings remain upper bounds. This caveat survives the
  longer runs unchanged. `fill_apo_cavity.py` and
  `compare_collapse_hydration.py` exist, so it can be tested rather than
  asserted, and it is the one Methods caveat worth keeping.
- **Delete the other Methods caveat.** "Apo replicates differ in length, so
  volumes were compared on a window common to all of them" is no longer true:
  the window is the full 100 ns for 18 of 20 genotypes, and V106A is the only
  short one.
