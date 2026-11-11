# Reproduction validation — 10 November 2026

- All 15 frozen submission files match their SHA-256 reference checksums.
- All 11 supported numerical rebuild targets complete in the isolated local `nnrti-prep` environment. The paper collection contains all 22 submitted figures and tables under their final publication labels, plus supported rebuilt candidates.
- All 21 tests pass, including native spreadsheet recalculation and independent verification of all 753 structural statistical comparisons against source replicate data. The structural workbook retains 10,732 live formulas and full-precision inputs.
- Preservation tests verify binary archive round trips, checksum rejection, source-change detection, retained prior MD copies and compatibility paths after workspace consolidation.
- The local legacy pmx compatibility stack reproduces all 114 saved phase/replicate BAR estimates from integrated work within reported rounding.

`reproduction-environment.json` records the tested local package versions. Each build writes its detailed input/output hashes and environment record to `generated/reproduction.json`, with verification logs under `generated/validation/`. Curated artwork and workbooks remain available through exact export; [ARTIFACTS.md](ARTIFACTS.md) documents their sources and manual composition steps.
