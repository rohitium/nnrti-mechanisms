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

## 5. HIV-specific MD / free-energy resistance studies (by target)

### 5a. Reverse transcriptase / NNRTIs — the closest precedents to this work
The Jorgensen-lineage RT/NNRTI FEP papers computed ΔΔG for the *same* NNRTI-pocket
residues in our panel (K103N, V106A, Y181C, L100I) with nevirapine/efavirenz, using
the same core logic (ΔΔG_bind vs observed resistance). These are our nearest
literature comparison and belong in the Discussion.

- **Energetic effects for observed and unobserved HIV-1 RT mutations of residues
  L100, V106, and Y181 in the presence of nevirapine and efavirenz.** MC/FEP.
  https://www.wikidata.org/wiki/Q33982889
  Key parallel: the clinically-*observed* variant had the more positive ΔΔG than
  unobserved codon alternatives — i.e. ΔΔG_bind tracked which mutation survives in
  patients. Direct precedent for our fold-change hypothesis test.
- **HIV-1 RT variants Y181C, V106A, L100I, K103N with NNRTIs — MC + linear-response
  (2004).** https://pubmed.ncbi.nlm.nih.gov/15553926/
  Same four residues we study; linear-response binding estimates per mutation.
- **Structural and Energetic Analyses of the K103N Mutation on Efavirenz Analogues.**
  *J. Med. Chem.* 2004. doi:10.1021/jm0303507 — https://doi.org/10.1021/jm0303507
  and **Activity predictions for efavirenz analogues with K103N (2003)**,
  https://pubmed.ncbi.nlm.nih.gov/12951121/ — MC/FEP on K103N (a charge-changing
  mutation, like our K103N/G190E legs); relevant to how they handled it.
- **Energetics of Mutation-Induced Changes in Potency of Lersivirine against HIV-1
  RT.** MC/FEP on a later NNRTI.
  https://www.academia.edu/85831907/ — extends the same framework to a
  next-generation NNRTI (closer in chemotype to doravirine).
- **MD of HIV-1 RT–DNA–nevirapine complexes explains NNRTI inhibition and resistance
  by connection-domain mutations (2013).** https://pubmed.ncbi.nlm.nih.gov/24174331/
  Resistance from mutations *outside* the pocket — a mechanism our pocket-only ΔΔG
  cannot capture (useful for the limitations section).
- **Thermodynamics of HIV-1 RT in Action.** *JACS* 2013.
  doi:10.1021/ja4018418 — https://doi.org/10.1021/ja4018418
- **Structural basis for drug resistance mechanisms for NNRTIs of HIV RT (2008).**
  https://pubmed.ncbi.nlm.nih.gov/18313784/ — structural context for pocket mutations.

### 5b. Protease — the DRV/decomposition lineage (our seed reference's neighbours)
- **Decomposing the energetic impact of drug-resistant mutations in HIV-1 protease
  on binding DRV.** https://pmc.ncbi.nlm.nih.gov/articles/PMC2882104/  (the seed paper)
- **Computational Mutation Scanning and Drug-Resistance Mechanisms of HIV-1 Protease
  Inhibitors.** https://pmc.ncbi.nlm.nih.gov/articles/PMC2916083/
  FEP absolute binding of APV and DRV consistent with experiment; decomposition shows
  mutations distort the active site so lost binding free energy is *not* confined to
  the mutated residues — a caution for per-residue interpretation.
- **Structural/dynamic/thermodynamic basis of DRV resistance in a heavily mutated
  protease (2022).** https://pmc.ncbi.nlm.nih.gov/articles/PMC9420863/
- **M46I-induced saquinavir resistance — MD + binding-energy (2022).**
  https://pmc.ncbi.nlm.nih.gov/articles/PMC9031992/
- **Susceptibility of HIV-1 protease variants to DRV and KNI-1657 — multiple MD +
  MM-PBSA/SIE.** *Langmuir* 2021. doi:10.1021/acs.langmuir.1c02348 —
  https://doi.org/10.1021/acs.langmuir.1c02348
- **FLAP+ / ACT variants vs amprenavir & darunavir — MM-PBSA (2015).**
  *Sci. Rep.* 5, 10517. https://www.nature.com/articles/srep10517

### 5c. Integrase / INSTIs
- **Mechanisms of HIV-1 integrase resistance to dolutegravir (2023).**
  *Sci. Adv.* doi:10.1126/sciadv.adg5953 — https://doi.org/10.1126/sciadv.adg5953
  All-atom free-energy simulations explain G140A/Q148K resistance via Mg²⁺
  polarization and weakened chelation — recent, mechanistically rigorous.
- **HIV-1C integrase E92Q/G140S/Y143R vs dolutegravir — MD binding free energy
  (2019).** *PLOS ONE.* https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0223464
- **Cross-resistance to INSTIs — MD + residue interaction network.** *J. Chem. Inf.
  Model.* https://pubs.acs.org/doi/abs/10.1021/ci300541c

### 5d. Capsid / lenacapavir
- **MD Free-Energy Simulations Reveal the Mechanism of M66I Capsid Resistance (2021).**
  *Viruses* 13(5), 920. doi:10.3390/v13050920 — https://doi.org/10.3390/v13050920
  Striking mechanistic parallel to our F227C reasoning: resistance comes from the
  *free-energy cost of side-chain reorganization* (I66 steric clash), NOT reduced
  protein–ligand interaction — exactly the kind of endpoint/reorganization effect that
  a bound-state-only ΔΔG can miss.
