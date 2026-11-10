# JCIM revision upload files

- **a. Point-by-point response:** Response-to-Reviewers.docx
- **b. Clean manuscript:** Manuscript-Clean.docx
- **c. Marked manuscript:** Manuscript-Marked.docx — upload as Supporting Information for Review. Changes relative to the submitted manuscript are shown in red.
- **d. Clean Supporting Information for Publication:** Supporting-Information-Clean.pdf. An editable DOCX is also provided.
- Upload Supplementary Tables 1–8 as the accompanying publication XLSX files. Their numbers are unchanged.
- Supporting-Information-Cover-Page-Clean.docx is the optional clean file inventory.

The clean and marked manuscripts have identical text, equations, tables, citation instructions and figure content. The clean SI matches the corrected marked SI. Revision colors, highlighting and tracked changes are absent from the clean copies. Embedded raster artwork was recovered from EMF wrappers without redrawing the figures so all figures render in PDF.

## Numerical audit

8,464 cell/statistical consistency checks passed, including independent reconstruction of all 753 Supplementary Table 4 statistics rows from full-precision replicate summaries, main Tables 1–3 versus SI, MM/GBSA and alchemical leg aggregation, and both main heatmap value/error matrices. All numeric source cells in Supplementary Tables 7–8 were also compared. Energy–phenotype regressions and classification counts were independently confirmed. Y318F phenotype provenance was verified against Table 2 of Asante-Appiah et al. (2021), DOI 10.1128/AAC.01216-21, because its individual isolate values are not published in the workbook.

Corrections include the G190S and Y318F induced-fit p-values (both 0.02 at two decimal places), two V106 comparison confidence intervals, one contact interval rounding endpoint, Y318F's non-zero Val106-distance SEM, protein-level caption counts, and the finite-size residual/uncertainty range. Workbook confidence limits now use full-precision Student t quantiles; p-values in the workbook were already correct. The Table 4 explanatory notes and generator were synchronized. Phenotype values originally reported at whole-number precision were retained; these are rounded displays, and unpublished Y318F raw values were not invented.

The audit verifies consistency of the supplied calculations and reporting; it does not rerun the molecular-dynamics simulations or establish their scientific convergence.

TOC graphic: TOC-Graphic.pdf and TOC-Graphic.tif (600 dpi). The identical original artwork is already on the final page of both manuscript copies, labeled TOC Graphic, without a caption. Sources: ops/figures/toc_graphic.pml and src/nnrti/cli/plot_toc_graphic.py.
