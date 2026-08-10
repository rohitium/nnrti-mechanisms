# Related work: FEP / alchemical free energy for drug-resistance mechanisms

Curated references for the DOR–RT non-equilibrium FEP manuscript, grouped by
relevance to our protocol (pmx hybrid topology + GROMACS non-equilibrium
switching + BAR, computing ΔΔG_bind per resistance mutation). Compiled 2026-08-08.

**Read these three first:** Gapsys 2020 (our exact method), Hauser 2018 (the
canonical resistance-by-FEP precedent), Serra 2025 (dissipation — speaks straight
to our charge-leg non-convergence and the F227C 500 ps rerun).

---

## 1. Protocol lineage — pmx / non-equilibrium alchemy

- **Gapsys V, Pérez-Benito L, Aldeghi M, Seeliger D, van Vlijmen H, Tresadern G,
  de Groot BL (2020).** "Large scale relative protein–ligand binding affinities
  using non-equilibrium alchemy." *Chem. Sci.* 11(4), 1140–1152.
  doi:10.1039/c9sc03754c — https://doi.org/10.1039/c9sc03754c
  Canonical validation of exactly what we run (pmx + GROMACS NEQ switching + BAR):
  482 perturbations across 13 datasets, AUE 3.64 kJ/mol, on par with FEP+. Cite as
  the provenance of the protocol in Methods.

- **Serra E, Ghidini A, Decherchi S, Cavalli A (2025).** "Nonequilibrium Binding
  Free Energy Simulations: Minimizing Dissipation." *J. Chem. Theory Comput.*
  21(4), 2079–2094. doi:10.1021/acs.jctc.4c01453 —
  https://doi.org/10.1021/acs.jctc.4c01453
  Dissipation is THE precision-limiting quantity in NEQ. Directly supports our
  charge-leg diagnosis (K103N/G190E non-convergence) and the F227C 500 ps rerun
  (slower switching → less dissipation → better overlap).

- **de Groot lab et al. (2025).** "On Free Energy Calculations in Drug Discovery."
  *Acc. Chem. Res.* doi:10.1021/acs.accounts.5c00465 —
  https://doi.org/10.1021/acs.accounts.5c00465
  Current state-of-the-field perspective; good single citation for the intro.

## 2. Flagship "FEP predicts resistance" applications

- **Hauser K, Negron C, Albanese SK, Ray S, Steinbrecher T, Abel R, Chodera JD,
  Wang L (2018).** "Predicting resistance of clinical Abl mutations to targeted
  kinase inhibitors using alchemical free-energy calculations." *Commun. Biol.*
  1, 70. doi:10.1038/s42003-018-0075-x —
  https://doi.org/10.1038/s42003-018-0075-x
  The reference resistance-by-FEP paper: 144 clinical Abl mutations, 8 inhibitors,
  RMSE 1.3 kcal/mol, ~88–93% resistant/susceptible classification. Our strongest
  precedent that ΔΔG_bind can rank resistance — cite where we frame the
  fold-change comparison as a hypothesis test.

- **Woods CJ, Malaisree M, Long B, McIntosh-Smith S, Mulholland AJ (2012).**
  "Quantitative Predictions of Binding Free Energy Changes in Drug-Resistant
  Influenza Neuraminidase." *PLoS Comput. Biol.* 8(8), e1002665.
  doi:10.1371/journal.pcbi.1002665 —
  https://doi.org/10.1371/journal.pcbi.1002665
  Two points we rely on: (a) resistance mutations give *subtle* 1–3 kcal/mol
  shifts — our error regime; (b) different mutations resist via different
  conformational routes, so a single mechanism is unlikely. Cite when defending
  divergence-as-finding / resistance-is-multi-mechanistic.

- **Predicting Resistance to Small Molecule Kinase Inhibitors (2025).**
  *J. Chem. Inf. Model.* Recent Abl-style extension; most up-to-date methodology
  comparison for resistance ranking.

## 3. HIV-specific free-energy studies

- **HIV-1 protease / darunavir ΔΔG decomposition.** "Decomposing the energetic
  impact of drug-resistant mutations in HIV-1 protease on binding DRV."
  https://pmc.ncbi.nlm.nih.gov/articles/PMC2882104/
  Closest in-virus precedent: per-mutation ΔΔG decomposition on the other major
  HIV target.

- **HIV-1 capsid M66I (2021).** "MD Free Energy Simulations Reveal the Mechanism
  for the Antiviral Resistance of the M66I HIV-1 Capsid Mutation."
  https://pmc.ncbi.nlm.nih.gov/articles/PMC8156065/
  Free-energy simulation resolving a resistance mechanism in HIV beyond RT/protease.

- **DOR clinical/structural context (not FEP, for the intro).** Doravirine
  susceptibility: >30× better than efavirenz vs K103N; Y188L drives complete
  resistance; Y181C partial. Real-world isolates:
  https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8597775/ ; IAS-USA 2025 mutations
  update: https://www.iasusa.org/wp-content/uploads/2025/03/33-2-mutations.pdf

## 4. Other viral targets — recent predict-then-validate templates

- **Discovery of Nirmatrelvir Resistance Mutations in SARS-CoV-2 3CLpro: A
  Computational–Experimental Approach (2023).** *J. Chem. Inf. Model.* 63.
  doi:10.1021/acs.jcim.3c01269 — https://doi.org/10.1021/acs.jcim.3c01269
  FEP alanine scanning predicted resistance mutations later confirmed by IC50
  (8–72× shifts). Clean template for the predict-then-validate narrative.

- **Nirmatrelvir P132H / P132H-A173V (2024).** *J. Chem. Inf. Model.*
  doi:10.1021/acs.jcim.4c00334 — https://doi.org/10.1021/acs.jcim.4c00334
  Current bar for combining alchemical free energy with mechanism dissection on a
  protease; useful if a reviewer pushes for more mechanism.
