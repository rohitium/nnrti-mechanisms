# Manuscript Feedback Notes

This file records collaborator feedback on the manuscript draft as objectively as possible. It combines:

- Comments and highlighted text extracted from [DorDRM-MD-04-13-26.docx](/Users/rohitpro/Career/00_Github/nnrti-mechanisms/manuscript/DorDRM-MD-04-13-26.docx)
- Notes from the follow-up call

## Source Context

- Collaborator comments in the Word document are attributed to Robert W. Shafer and are dated April 14-15, 2026.
- The Word file contains 12 comment balloons and 25 yellow-highlighted passages.
- Some highlighted passages have explicit written comments attached; others are highlighted without a separate written note.

## Call Notes: Overall Direction

- The overall flow of the manuscript needs work, especially the transitions between binding energy calculations, residue contact analysis, and logistic regression training/testing.
- It is currently unclear how one section motivates or leads into the next.
- It is unclear whether the present section ordering is the right overall flow for the manuscript.
- The Introduction has a completeness problem.
- The manuscript should cite all recent real-world evidence and clinical trial data regarding doravirine.
- `Rhee et al.` should remain in the references, but should not be used as a primary data source.
- For Figure 1, the fold-change resistance data for each genotype should be individually traced to and cited from the original publication source.
- The binding energy section should more explicitly describe how the calculations are performed.
- The binding energy section should explain why the current calculations do not represent absolute binding free energies.
- The binding energy section should explain what would be required to calculate the RT-DOR binding energy correctly from equilibrium simulations.
- In the structural mechanism section, residues such as Pro225 for Y188L and Ser105 for V106A currently appear arbitrary to the reader.
- The manuscript should explain what motivated examination of those residues.
- The manuscript should clarify whether residue/contact selection followed a systematic process rather than post hoc cherry-picking.
- Structural descriptions such as `"P225H tilting Tyr181 toward rim of the pocket rather than DOR"` or `"V106A causes shift in DOR pose"` are currently too abstract.
- Figures and/or quantitative analyses should show these structural changes more clearly.
- Figure 3 should include occupancies of all residues that ever contact DOR.
- The occupancy threshold of `mean_occupancy = 0.5` should be stated explicitly as a tractability choice.
- The logistic regression section should show the full model details so the meaning of `"predicted probability of resistance"` is explicit.
- The logistic regression section should explain how class labels and the model threshold of `0.5` affect the reported metrics.
- The term `"test set"` should be changed to `"Limited data set"`.
- The rationale for renaming `"test set"` is that these data are not necessarily fully accurate and are being interrogated by the simulations/modeling rather than serving as a clean benchmark.
- The manuscript should emphasize that the model uses only 3 features.
- The manuscript should emphasize that the model predicts low probability of resistance for WT as an unseen case.

## Extracted Word Comments

### Citation and Source Coverage

- Anchor: `"Mejias-Trueba et al., 2024"`
- Comment: `"Add others - e.g.: https://www.sciencedirect.com/science/article/pii/S2352301824001504"`

- Anchor: `"Monogram PhenoSense platform"`
- Comment: `"And the Merck reporter gene assay."`

- Anchor: `"Stanford HIVDB to obtain phenotypic susceptibility data"`
- Comment: `"Consider using PhenoSense and Merck"`

- Anchor: `"Rhee et al., 2023"`
- Comment: `"Dont cite us as a primary source."`

### Results Framing and Data Provenance

- Anchor: `"In the Monogram-based analysis of 4,070 clinical isolates, DOR retained activity against more isolates than EFV, nevirapine (NVP), rilpivirine (RPV), or etravirine (ETR), depending on the biological cutoff applied (Asante-Appiah et al., 2021)"`
- Comment: `"Indicate that this analysis was based on mutations containing a single mutation."`

### Binding Energy Terminology and Interpretation

- Anchor: `"in Electrostatic Binding Energy"`
- Comment: `"Is this the same as delta-deltaG"`

- Anchor: `"loss of this interaction in case of Y188L leads to not only a significant reduction in binding energy (Supplementary Figure 1) but also alteration in pocket geometry that transiently shifts Pro225 away from DOR "`
- Comment: `"Why isn't the reduction in binding energy enough?"`

### Structural Interpretation and Clarity

- Anchor: `"is tilted toward the rim of the pocket, rather than DOR"`
- Comment: `"Unclear"`

- Anchor: `"Pro225 is located near the entrance of the NNIBP in the primer grip (β12–β13 hairpin) region, which means that this remodeling likely increases the risk of DOR escaping the NNIBP in the Y188L genotype."`
- Comment: `"But Y188L and P225H dont occur together"`

### Logistic Regression Framing

- Anchor: `"accuracy of 89.5% (100% on train and 75% on test set) and ROC AUC of 83.3% (100% on train and 41.7% on test set"`
- Comment: `"Implies a threshold"`

- Anchor: `"However, this sample size is quite small and the test ROC AUC metric of 41.7% strongly suggests this model does not yet generalize well to unseen data."`
- Comment: `"I thought the AUC was much higher."`

### Scope and Limitations

- Anchor: `"Fourth, detailed structural mechanisms of resistance to two canonical DOR DRMs, Y318F and A98G+F227C, were not elucidated from our work"`
- Comment: `"I thought they were not included"`

## Highlighted Passages Without Separate Written Comment

These passages are highlighted in yellow in the Word document but do not have a separate written comment attached.

### Abstract and High-Level Claims

- `"demonstrates weak but statistically significant correlation with the available phenotypic DOR susceptibility changes in vitro"`
- `"producing probabilities that correlate well with available in vitro DOR susceptibility data."`

### Results, Binding Energy, and Interpretation

- `"M230L"`
- `"only"`
- `"∆∆G"`
- `"Moreover, most components of the binding energy computed here do not correlate with the phenotypic susceptibility data either (Supplementary Figure 1), except the electrostatics component, which demonstrates weak but statistically significant correlation (Figure 2)."`
- `"an accurate estimate of binding energy would require sampling equilibrium populations of bound and unbound complexes, lack of correlation between the binding energy and observed phenotypic susceptibility is not surprising"`

### Contact Metrics and Occupancy Analysis

- `"heavy-atom distance"`
- `"for all residues"`
- `"which were in contact with DOR at least 50% of the time, i.e. Mean Occupancy ≥ 0.5, in any simulation"`

### Y188L Mechanism

- `"Y188L preferentially alters binding pocket geometry, causing Pro225 to shift away from DOR."`
- `"showing Pro225 almost twice the distance away from DOR in case of Y188L"`

### V106A and Related Mechanisms

- `"V106A preferentially alters DOR pose, causing Ser105 to shift closer to DOR"`
- `"V106A has a significant impact on the interaction with the central pyridone core in DOR"`
- `"Similar structural shift toward the pocket entrance"`
- `"Leu"`
- `"Ser105"`
- `"106x vs 105x"`
- `"Val179"`

### Logistic Regression and Classification Language

- `"DOR-Ser105 Distance, DOR-Tyr188 Distance, and the Root Mean Square Deviation (RMSD) of the heavy atoms in DOR from the baseline conformation in the crystal structure"`
- `"“high” resistance"`

## Consolidated Themes

- Literature coverage and source attribution are incomplete.
- Figure 1 requires genotype-by-genotype source validation from original publications.
- The manuscript needs clearer definitions around electrostatic energy, `∆∆G`, and what is or is not being estimated by the MD-based calculations.
- The manuscript needs a more explicit and systematic explanation of how residue contacts were selected and interpreted.
- Structural mechanism claims need clearer visual or quantitative support.
- Occupancy analysis should be presented more comprehensively and with explicit thresholding rationale.
- The logistic regression section needs fuller methodological transparency and more careful language around labels, thresholds, probabilities, and evaluation.
- Several current claims may read as stronger or more definitive than the collaborator is comfortable with.
