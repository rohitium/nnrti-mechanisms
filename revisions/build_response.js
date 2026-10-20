const fs = require('fs');
const { Document, Packer, Paragraph, TextRun, HeadingLevel, AlignmentType } = require('docx');
const B = JSON.parse(fs.readFileSync('/private/tmp/claude-502/-Users-rohitpro-Career-00-Github-nnrti-mechanisms/73e5c257-8937-4941-835d-14bb0d9b9bec/scratchpad/reviewer_blocks.json', 'utf8'));

const RED = "FF0000";
const kids = [];

const title = (t) => kids.push(new Paragraph({ heading: HeadingLevel.HEADING_1, spacing:{after:200},
  children:[new TextRun({ text:t, bold:true })] }));
const h2 = (t) => kids.push(new Paragraph({ heading: HeadingLevel.HEADING_2, spacing:{before:320, after:160},
  children:[new TextRun({ text:t, bold:true })] }));
// reviewer text, verbatim, black and italic
const quote = (i) => kids.push(new Paragraph({ spacing:{after:140}, indent:{left:420},
  children:[new TextRun({ text:B[i], italics:true })] }));
// our response, red
const resp = (...ps) => ps.forEach((t, n) => kids.push(new Paragraph({ spacing:{after:160},
  children:[
    ...(n===0 ? [new TextRun({ text:"Response: ", bold:true, color:RED })] : []),
    new TextRun({ text:t, color:RED })
  ]})));
const pending = (t) => kids.push(new Paragraph({ spacing:{after:200}, indent:{left:420},
  children:[new TextRun({ text:"[PENDING — "+t+"]", bold:true, color:RED })] }));

title("Response to Reviewers");
kids.push(new Paragraph({ spacing:{after:120}, children:[new TextRun({
  text:"Manuscript ci-2026-030865 — “Molecular Simulation of Doravirine Resistance in HIV-1 Reverse Transcriptase”", bold:true })]}));
kids.push(new Paragraph({ spacing:{after:260}, children:[new TextRun({
  text:"Rohit Satija, Kaiming Tao, Robert W. Shafer" })]}));

resp("We thank the editor and the three reviewers for a careful and constructive set of reports. The reviewers converge on one substantive point — that the sampled DOR-bound ensemble is treated as representative without that assumption being tested — and we have taken it seriously, both in new analysis and in how the claims are now framed. Reviewer comments are reproduced below in black; our responses are in red. All changes to the manuscript are marked in red in the revised file.");

// ---------------- Reviewer 1 ----------------
h2("Reviewer 1");
[0,1,2,3,4].forEach(quote);
resp("We appreciate this framing, which identifies the central limitation more precisely than we had. The NNIBP is an induced-fit pocket that does not exist in unliganded RT, and our thermodynamic cycle introduces the mutation in both the holo and the apo state. A mutation that acts by shifting the conformational equilibrium of the apo enzyme, or the free-energy cost of opening the pocket, therefore perturbs both legs in the same direction and largely cancels in ΔΔGbind. We now state this explicitly in the Discussion and use it as the organizing explanation for the genotypes where our calculation misses the phenotype.");

h2("Reviewer 1, comment 1 — weak correlation with phenotype");
[5,6,7].forEach(quote);
resp("We agree and have tempered the claim accordingly. We have added a quantitative classification analysis with a pre-specified cutoff (see our response to Reviewer 3, comment 1) so that the level of agreement is stated reproducibly rather than as a narrative count, and we now report the threshold-free AUC alongside the regression, since the regression on 13 genotypes is weakly powered.",
"On the substantive point, we agree that the weak correlation is unlikely to be purely a free-energy-calculation error. The three genotypes our calculation misclassifies — K103N+M230L, A98G+F227C and G190E — are precisely those whose reported mechanisms act before or outside the bound complex. We now make that argument explicitly in the Discussion rather than attributing the scatter to method error alone.");

h2("Reviewer 1, comment 2 — convergence within a basin");
[8,9,10].forEach(quote);
resp("We accept this distinction and have adopted the reviewer's language in the Discussion: stable RMSD values, overlapping work distributions and small replicate errors establish convergence within the sampled basin, not that the basin is the relevant one. We also now state plainly that our three replicates were seeded from the same crystallographic pose with only a 0.1 Å coordinate perturbation, so they are not independent in the sense the reviewer means.",
"On the apo simulations, we note that our existing apo data is more extensive than the submitted manuscript conveyed. We hold three 100 ns apo replicates for A98G+F227C, F227C and G190E, and 37–75 ns for K103N+M230L, which covers all three genotypes where ΔΔGbind misses the phenotype as well as the M230-containing genotype the reviewer names. We are completing a WT-referenced comparison of the apo and holo pocket ensembles for these systems.");
pending("apo-vs-holo ensemble comparison (pocket volume, RMSF, 225–230 and 100–110 loop conformations) and the independent-starting-pose test. Insert results and state which simulations were added, if any.");

h2("Reviewer 1, comment 3 — the mismatches suggest different mechanisms");
[11,12,13].forEach(quote);
resp("We agree that these are not random errors, and the revised Discussion treats them as mechanistically informative. For G190E and A98G+F227C, the interface distortion is real in the holo state (ΔΔETotal = 1.94 ± 0.40 and 2.47 ± 0.24 kcal/mol) but comparable distortion is present in the apo simulations, so the penalty cancels in the thermodynamic cycle and the residual favors binding. That cancellation is the mechanism, not an artifact to be explained away, and it is exactly what the reviewer anticipates when asking whether the remodeling occurs in both apo and holo RT.",
"For K103N+M230L the DOR interface is essentially unperturbed in our trajectories — chlorocyanophenyl burial is the second highest in the panel and the Lys103 anchor is intact — while Met230 lies in the p66 primer grip and contacts the primer terminus. Our binary complex contains no template/primer, so this mechanism is structurally absent from the simulation rather than poorly sampled by it.");

h2("Reviewer 1, comment 4 — clearer summary of the structural analysis");
[14,15,16].forEach(quote);
resp("We agree that the structural results are hard to navigate as three dense tables. We are adding a figure mapping all 18 genotypes onto the RT–DOR structure relative to DOR and the NNIBP, and a WT-referenced heatmap of the Table 3 metrics, and are extending the analysis to the protein-level observables the reviewer lists — RMSF, secondary-structure changes and pocket conformational distributions — so that direct changes in DOR contacts can be separated from broader remodeling of the pocket.");
pending("new mutation-map figure, WT-referenced heatmap, and RMSF / secondary-structure / pocket-distribution analyses. Insert figure numbers and findings.");

// ---------------- Reviewer 2 ----------------
h2("Reviewer 2");
[17,18,19,20].forEach(quote);
resp("We thank the reviewer for the accurate summary and for the methodological points below, several of which identified real gaps in what we reported.");

h2("Reviewer 2, comment 1 — convergence diagnostics");
quote(21);
resp("Block averaging over the production trajectories was already part of our analysis but was not reported. We are now including it together with autocorrelation times and the resulting effective sample sizes, applied to the structural observables and to the MM/GBSA and FEP estimates, and we discuss in the Discussion how the limited sampling may bias the results — specifically that our replicates establish convergence within one conformational basin rather than across the accessible landscape.");
pending("autocorrelation / effective-sample-size analysis. Insert the table and the resulting statement about sampling adequacy.");

h2("Reviewer 2, comment 2 — the binary RT–DOR system");
quote(22);
resp("We agree, and the revised Discussion states explicitly which mechanisms this setup forecloses. The binary complex contains no template/primer and no incoming dNTP, so resistance acting through the primer grip (M230L) is structurally inaccessible to our calculations, as is any mechanism operating through the catalytically productive ternary complex. We also note that resistance can act through the conformational equilibrium of the unliganded pocket on timescales far longer than we sample, and through effects on processivity and replicative capacity, none of which enter ΔΔGbind. We argue that this is a leading explanation for the K103N+M230L false negative rather than a generic caveat.");

h2("Reviewer 2, comment 3 — effect sizes and confidence intervals");
quote(23);
resp("We accept that several mechanistic statements rested on differences presented only as means with replicate SEMs. We are adding bootstrap confidence intervals on the effect sizes for the specific comparisons the reviewer names, V106A versus V106I/M and G190A versus G190E, and are supporting the ligand-displacement claims with a torsion-level analysis of the four rotatable DOR dihedrals in addition to the per-moiety contact counts already reported.");
pending("bootstrap confidence intervals and the DOR torsion analysis. Insert effect sizes with CIs and the torsion distributions.");

h2("Reviewer 2, comment 4 — charge-changing legs");
quote(24);
resp("This prompted us to verify rather than assume, and we can now state the result explicitly in the Methods. Protonation states are identical between the wild type and every mutant by construction: all systems inherit their states from the same prepared wild-type structure, in which hydrogens were built once and carried unchanged through mutation modeling. We confirmed this directly — all thirteen histidines are neutral with the proton on Nδ1, identically across the wild type and all 18 genotypes — so no titratable residue in the NNIBP differs between the end states of any leg. The glutamate introduced by G190E was modeled deprotonated and singly charged, which is what makes that leg charge-changing.",
"The Methods now also describe the hybrid residue topologies (pmx single-topology construction, with atoms unique to one end state present as non-interacting dummies at the other, bonded parameters interpolated linearly in λ, dummies relaxed during end-state minimization) and give the soft-core parameters as run: the linearized Gapsys form with α = 0.3, σ = 0.25 nm, soft-core power 1, applied to both Lennard-Jones and Coulomb interactions.",
"On the finite-size correction we have chosen to be explicit rather than to imply completeness. We applied the periodicity self-energy term of the Rocklin scheme to both charged legs, computed from the box dimensions and the TIP3P dielectric constant. The higher-order terms — net-charge undersolvation and discrete-solvent — were not applied in the submitted work, and we now state this and report their estimated magnitude so the residual uncertainty on these two legs is visible to the reader.");
pending("magnitude estimate for the omitted higher-order Rocklin terms, or their implementation. Insert the numbers for the K103N and G190E legs.");

h2("Reviewer 2, comment 5 — Gasteiger charges");
quote(25);
resp("We accept this and have done both of the things the reviewer offers as alternatives. The Discussion now explains why the choice affects our two calculations unequally: MM/GBSA evaluates an absolute interaction energy, so ΔEelec and ΔEGB depend directly on the ligand charge set, whereas the alchemical calculation perturbs the protein and not the ligand — DOR is neutral and parameterized identically at both end states of every leg, so the ligand contribution largely cancels in the thermodynamic cycle. We are also recomputing MM/GBSA with AM1-BCC charges for a subset of genotypes to quantify the effect directly rather than argue it.");
pending("AM1-BCC MM/GBSA recomputation for a genotype subset. Insert the comparison.");

// ---------------- Reviewer 3 ----------------
h2("Reviewer 3");
[26,27,28].forEach(quote);
resp("We thank the reviewer, and we agree that the central quantitative claim needed either firmer statistical support or more careful wording. We have taken both routes: the agreement with phenotype is now stated against a pre-specified threshold with a full confusion matrix and a threshold-free AUC, and the claim itself is tempered.");

h2("Reviewer 3, comment 1 — define “correctly predicted”");
quote(29);
resp("Agreed, and now stated. A genotype is counted as predicted to reduce binding when ΔΔGbind exceeds +0.5 kcal/mol, approximately the thermal energy kBT at 300 K, and as predicted to have little effect otherwise; the same cutoff is applied to ΔΔETotal. The cutoff was fixed a priori from the thermal energy rather than fitted to the panel. Against the 13 genotypes with an established phenotype, ΔΔGbind gives 7 true positives, 3 true negatives, 1 false positive (V106I) and 2 false negatives (A98G+F227C and K103N+M230L): sensitivity 0.78, specificity 0.75, Matthews correlation coefficient 0.50. Because any single cutoff is arbitrary for an observable that was never calibrated to a clinical threshold, we also report the area under the ROC curve, 0.78 (exact permutation p = 0.074). This appears in the Abstract, in the Statistical analysis section and in the Results.");

h2("Reviewer 3, comment 2 — typographical issues");
quote(30);
resp("We have proofread the manuscript in full and corrected the errors found. On the three specific items: we were unable to locate the strings “NRNTI-resistance mutations” or an “APPROVAL” section header in our submitted manuscript file or Supporting Information, and we suspect these arose in the PDF generated during submission rather than in the source document. We have checked the section headings, and the author-information block contains only Author Contributions, Funding Sources and Notes. The manuscript has additionally been reviewed in full for language and usage.");

h2("Reviewer 3, comment 3 — MM/GBSA minimization protocol");
[31,32,33].forEach(quote);
resp("(a) The justification is now given in the Methods, with citations. Minimizing each snapshot before the energy evaluation is established practice for end-point methods rather than a departure from it: minimization-based MM-PBSA was validated for drug discovery by Kuhn et al. (J. Med. Chem. 2005, 48, 4040), Rastelli et al. (J. Comput. Chem. 2010, 31, 797) showed that minimization-based MM-PBSA and MM-GBSA reproduce the accuracy of ensemble averaging at a fraction of the cost, and minimizing from a series of MD snapshots \u2014 which is what we do \u2014 was shown by Weis et al. (J. Med. Chem. 2006, 49, 6596) and in the review we already cite (Genheden and Ryde, ref 34) to give results similar to single-point ensemble averages. That review states directly that minimized structures often give as good or better results as MD simulations.",
"The purpose of the minimization in our case is to relieve steric strain that single-point evaluation of frames sampled under explicit solvent, periodic boundary conditions and constrained bond lengths would otherwise carry into the cutoff-free, unconstrained implicit-solvent Hamiltonian used for scoring, where it contributes Lennard-Jones energies far larger than the binding signal. We agree with the reviewer that this needed stating and defending rather than simply reporting, and the revised Methods now does so.",
"(b) We accept the reviewer\u2019s point that full minimization may reduce energetic noise and bias the components, and we have recomputed a subset of genotypes (WT, K103N, V106A, Y188L, G190E and V106I+F227C, three replicates each) under single-point evaluation, under hydrogen-only relaxation with heavy atoms restrained, and under the unrestrained full minimization used in the main text, reporting \u0394\u0394ETotal and each component side by side so the sensitivity is visible (Supplementary Table 5).");

pending("insert the Supplementary Table 5 numbers and one sentence stating whether the conclusions are unchanged.");

const doc = new Document({
  creator: "Rohit Satija",
  title: "Response to Reviewers — ci-2026-030865",
  sections: [{ properties:{ page:{ size:{ width:12240, height:15840 } } }, children: kids }]
});
Packer.toBuffer(doc).then(b => {
  fs.writeFileSync('/private/tmp/claude-502/-Users-rohitpro-Career-00-Github-nnrti-mechanisms/73e5c257-8937-4941-835d-14bb0d9b9bec/scratchpad/Response-to-Reviewers.docx', b);
  console.log('written,', kids.length, 'paragraphs');
});
