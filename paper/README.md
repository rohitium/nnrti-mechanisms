# Manuscript materials

The [November 10 submitted package](../revisions/submission-ready/) is the reference version. [ARTIFACTS.md](ARTIFACTS.md) maps its figures and tables to data and build commands.

- `sources/Figures.pptx`: latest editable figure deck; compare with final embedded artwork before editing because some panels were replaced directly in Word.
- `sources/figure1/`: structure artwork and supporting exports for Figure 1.
- `tables/Table-2-energetics.csv`: numerical reference for the energetics table.
- `submission-checksums.json`: SHA-256 manifest for the frozen submitted files.
- `artifacts.json`: machine-readable registry of all 22 submitted figures and tables and their supported rebuild targets.
- `generated/paper/README.md` (ignored): consolidated paper collection, created by `../workflows/05_manuscript_artifacts.sh`.
- `generated/reproduction.json`, `generated/validation/`: provenance hashes, environment versions and test results.
- Other `generated/` directories: plotting intermediates. Safe to recreate.

Earlier manuscripts and original submission files are under the single root `archive/`. Submitted citation fields remain in the Word documents; exporting figures or tables does not alter them.
