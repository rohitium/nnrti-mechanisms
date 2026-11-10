# Cleanup and reproduction validation — 10 November 2026

- SHA-256 checks pass for all 15 submitted Word/PDF/Excel/TOC files before and after rebuilding. Submitted content was not edited.
- All 11 supported artifact targets complete from an isolated copy containing only files intended for GitHub, without local trajectories or checkpoints. Exact exports include six main and seven supplementary figures, their final captions, two main tables and all seven SI workbooks.
- All 14 tests pass, including native LibreOffice recalculation of the generated structural workbook. Independent NumPy/SciPy calculations verify all 753 statistical comparisons against source replicate data, both for the submitted cached workbook and the regenerated workbook. The generator retains full precision and writes 10,732 formula cells.
- Rebuilt Figures 3–6 and Supplementary Figure 7 were visually compared with their submitted artwork. No scientific analysis protocol or submitted manuscript wording was changed during cleanup.
- All historical directories were consolidated beneath root `archive/`; moves are recorded there. Large local trajectory backups remain ignored. No raw simulation files were deleted.

The tested Python package versions are recorded in `reproduction-environment.json`. Figure 1, Figure 2 and curated SI workbooks are preserved through exact export; the artifact map documents their component sources and remaining manual steps. Raw-data backup completeness on Sherlock was not re-audited during this cleanup.
