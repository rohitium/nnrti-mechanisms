# Switch to a 4.0 Å contact cutoff — every number that changes

Applied 2026-08-31. `CONTACT_CUT` in `compute_mechanism_coordinates.py` and
`CUTOFF_NM` in `compute_dor_moiety_contacts.py` are both now **4.0 Å**.
All derived files and figures regenerated at stride 1 (mechanism coordinates)
and stride 2 (moiety contacts).

Distances, angles, pocket volumes and energies are unaffected — no contact
threshold enters them. **Only contact counts change.**

---

## 1. Chlorocyanophenyl ring burial — Supplementary Figure 3A

| system | 4.5 Å (old) | **4.0 Å (new)** | ratio to WT, 4.5 | **ratio to WT, 4.0** |
|---|---:|---:|---:|---:|
| **WT** | 45.9 ± 1.3 | **19.5 ± 0.5** | 1.00× | 1.00× |
| **Y188L** | 34.6 ± 0.9 | **13.7 ± 1.5** | 0.75× | **0.70×** |
| **Y181C** | 48.6 ± 1.4 | **24.0 ± 2.6** | 1.06× | **1.23×** |
| K103N | 46.2 ± 1.7 | 20.3 ± 0.6 | 1.01× | 1.04× |
| G190A | 45.0 ± 0.3 | 20.0 ± 1.2 | 0.98× | 1.02× |
| G190S | 42.6 ± 1.8 | 19.0 ± 0.4 | 0.93× | 0.97× |
| G190E | 46.6 ± 1.9 | 22.7 ± 2.0 | 1.02× | 1.16× |
| V106A | 43.0 ± 1.4 | 17.8 ± 0.7 | 0.94× | 0.91× |
| V106I | 44.9 ± 1.6 | 20.2 ± 0.3 | 0.98× | 1.04× |
| V106M | 45.4 ± 2.1 | 19.3 ± 0.9 | 0.99× | 0.99× |
| K103N+M230L | 48.6 ± 1.3 | 22.8 ± 2.4 | 1.06× | 1.17× |
| K103N+P225H | 44.4 ± 1.3 | 18.3 ± 0.6 | 0.97× | 0.94× |
| L100I+K103N | 44.5 ± 1.2 | 18.4 ± 0.5 | 0.97× | 0.94× |
| A98G+F227C | *not computed* | **20.4 ± 3.5** | — | 1.04× |
| V106I+F227C | *not computed* | **19.2 ± 2.3** | — | 0.99× |
| V106A+F227L | *not computed* | **17.2 ± 0.3** | — | 0.88× |
| V106A+L234I | *not computed* | **18.9 ± 0.4** | — | 0.97× |
| V106A+P225H | *not computed* | **17.9 ± 1.2** | — | 0.92× |
| Y318F | *not computed* | **19.2 ± 0.2** | — | 0.98× |
| F227C | *not computed* | **20.1 ± 0.3** | — | 1.03× |

**Bonus: the panel is now complete.** The regeneration covers all 20 systems;
previously A98G+F227C, V106I+F227C, Y318F and F227C had never been run through
`compute_mechanism_coordinates.py`, which was a gap flagged in
`discussion-expansion-2026-08-31.md`. That gap is closed.

Y188L remains the clear outlier — 0.70×, against 0.88× for the next lowest
(V106A+F227L).

---

## 2. Wording changes required

### (a) Y188L — the loss is now 30%, not 25%

> "…falling from 45.3 ± 2.4 in WT to 34.2 ± 1.1 in Y188L (Supplementary
> Figure 3A), a loss of roughly a quarter."

becomes

> …falling from **19.5 ± 0.5** in WT to **13.7 ± 1.5** in Y188L (Supplementary
> Figure 3A), a loss of **nearly a third**.

(This also folds in the earlier stale-number correction — the draft's
`45.3 ± 2.4 / 34.2 ± 1.1` was never reproducible from the committed data.)

### (b) Y181C — "nearly identical" no longer holds

At 4.5 Å, Y181C and WT were indistinguishable (48.6 vs 45.9, 1.06×). At 4.0 Å
Y181C is **23% higher** than WT (24.0 ± 2.6 vs 19.5 ± 0.5). The difference is
**1.7 σ — not statistically significant** (per-replicate Y181C 19.8 / 28.6 /
23.6 against WT 20.5 / 19.2 / 18.8), but "nearly identical" is no longer the
right description.

> "In simulations of the Y181C genotype, DOR packing in the NNIBP is nearly
> identical to WT: 48.6 ± 1.4 RT heavy atoms within 4.5 Å of the
> chlorocyanophenyl ring against 45.9 ± 1.3 in WT."

becomes

> In simulations of the Y181C genotype, packing around the chlorocyanophenyl
> ring is **not reduced** relative to WT — **24.0 ± 2.6 heavy-atom contacts
> within 4.0 Å against 19.5 ± 0.5 in WT**, a difference that is if anything
> favourable and does not reach significance across replicates.

The conclusion is unchanged and arguably strengthened: Y181C does not cost DOR
any packing.

### (c) "RT heavy atoms" → "heavy-atom contacts"

Unrelated to the cutoff, but it must be fixed in the same edit: these are counts
of atom **pairs**, not atoms (see `discussion-expansion-2026-08-31.md`, second
correction). Every occurrence of "RT heavy atoms within 4.0 Å of X" should read
"RT heavy-atom contacts with X (atom pairs within 4.0 Å)".

---

## 3. Whole-ligand and moiety contacts (V106A set)

At 4.0 Å with the **symmetric** moiety definition (every ring keeps its own
exocyclic substituents — CF₃ and carbonyl on the pyridinone, methyl and carbonyl
on the triazolinone, Cl and nitrile on the phenyl):

| genotype | whole ligand | chlorocyanophenyl | **pyridinone** | triazolinone |
|---|---:|---:|---:|---:|
| WT | 101.9 ± 2.9 | 29.8 ± 0.5 | **14.7 ± 1.4** | 35.4 ± 3.2 |
| V106A | 97.6 ± 2.2 | 28.6 ± 0.3 | **11.7 ± 1.8** | 35.9 ± 2.8 |
| V106A+F227L | 92.4 ± 1.7 | 26.6 ± 0.8 | **10.3 ± 0.7** | 32.8 ± 0.9 |
| V106A+L234I | 97.4 ± 2.2 | 28.7 ± 0.1 | **10.5 ± 0.8** | 34.2 ± 1.0 |
| V106A+P225H | 96.8 ± 1.9 | 27.5 ± 0.4 | **11.8 ± 0.1** | 33.2 ± 0.6 |

Percent change vs WT:

| genotype | whole | chlorocyanophenyl | **pyridinone** | triazolinone |
|---|---:|---:|---:|---:|
| V106A | −4.2% | −4.0% | **−20.4%** | +1.4% |
| V106A+F227L | −9.3% | −10.7% | **−29.9%** | −7.3% |
| V106A+L234I | −4.4% | −3.7% | **−28.6%** | −3.4% |
| V106A+P225H | −5.0% | −7.7% | **−19.7%** | −6.2% |

The headline result is **unchanged and slightly stronger**: the pyridinone loses
20–30% in every V106A genotype, four to seven times the fractional loss of any
other moiety. The whole-ligand number remains the misleading one (−4.2% to
−9.3%).

The `224 ± 1 → 212 ± 1` sentence becomes approximately **`102 ± 3 → 98 ± 2`**,
but it should be replaced wholesale by the moiety-resolved text — see
`discussion-expansion-2026-08-31.md` §Suggested replacement, with the
triazolinone "gain" clause removed (it does not survive at 4.0 Å either:
+1.4% / −7.3% / −3.4% / −6.2%).

### Partner side, 4.0 Å

Losing: Val106 (2.39 → 0, mutated), Tyr318 (5.07 → 3.43), Phe227 (3.57 → 2.69),
Lys101 (1.03 → 0.58), Leu100 (2.40 → 1.98).
Gaining: Ser105 (0.03 → 2.49), Ala106 (0 → 1.80), Lys104 (0.09 → 1.20).

Same directional story as at 4.5 Å.

---

## 4. Figures regenerated

| figure | file | change |
|---|---|---|
| Supp. Fig. 3A / mechanism panel | `results/analysis/mechanisms/plots/mechanism_panel.png` | y-axis rescaled; panel B (Ser105 distances) **unchanged**, being distances |
| Figure 1B | `results/plots/figure1B_dor_schematic.pdf` / `.png` | rebuilt: title and footer removed, RT residue structures drawn, connectors dropped except the Lys103 hydrogen bond, contact counts now 4.0 Å |

Figure 1B residue counts (whole-ligand contacts, 4.0 Å): Tyr188 20.3, Tyr318
9.5, Val106 9.5, Trp229 8.4, Lys103 6.7, Phe227 6.5, Tyr181 4.3, Gly190 3.0.

**Note on Tyr181 in Figure 1B:** at 4.0 Å it registers 4.3 whole-ligand contacts
but **0.1 with the chlorocyanophenyl ring and 0.0 with the other two rings** —
its contacts are entirely with the ether/methylene linker. The figure caption
should therefore say Tyr181 makes no contact *with the aromatic rings*, not that
it makes no contact at all. The earlier claim of "0.4 contacts with the whole
ligand" was per-moiety and should not be used unqualified.

---

## 5. Files regenerated

- `results/analysis/mechanisms/mechanism_coordinates.csv` (all 20 systems)
- `results/analysis/mechanisms/mechanism_summary.csv`
- `results/analysis/mechanisms/plots/mechanism_panel.png`
- `results/analysis/mechanisms/dor_moiety_contacts_{summary,per_replicate}.csv`
- `results/analysis/mechanisms/dor_residue_contact_delta.csv`
- `results/analysis/mechanisms/figure1B_contacts_4A.csv`
- `results/plots/figure1B_dor_schematic.{pdf,png}`

Not regenerated because unaffected: `dor_key_contacts_timeseries_all_mutations.csv`
(distances), `pocket_volume_profiles.csv`, `com_distance_profiles.csv`,
`y188_interplanar_angle_190series.csv` (angles), all binding-energy tables.
