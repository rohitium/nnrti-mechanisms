# Doravirine resistance in HIV-1 reverse transcriptase

Code, derived data and manuscript materials for the JCIM revision submitted **10 November 2026** (ci-2026-030865). The study compares WT RT with 18 resistance-associated genotypes using equilibrium MD, nonequilibrium FEP and MM/GBSA. F227C is an additional alchemical intermediate.

Start with the [submitted manuscript and SI](revisions/submission-ready/). The [artifact map](paper/ARTIFACTS.md) connects the current figure/table numbers to source data and commands. [REPRODUCE.md](REPRODUCE.md) explains setup, rebuilding and recovery of simulation data.

```bash
conda env create -f environment.yml
conda activate nnrti-prep
pip install -e .
python -m nnrti.cli.manuscript --check --export
# Rebuild plots and numerical tables from deposited derived data:
./workflows/05_manuscript_artifacts.sh
pytest -q
```

Exports and rebuilds go to ignored `paper/generated/`. The submitted Word, PDF and Excel files are protected by SHA-256 checksums and are never build destinations. Figure 1, Figure 2 and some SI tables include curated composition or experimental inputs; their exact submitted versions can be exported without recreating that manual work.

| Directory | Purpose |
|---|---|
| `revisions/submission-ready/` | Frozen November submission: clean/marked manuscript, response, SI and seven workbooks |
| `paper/` | Artifact map, submission checksums, editable figure sources and generated exports |
| `src/nnrti/` | Scientific code and analysis CLIs |
| `data/`, `manifests/` | Experimental inputs, structures, preparation inputs and run provenance |
| `results/analysis/` | Derived numerical data and plots; source filenames may retain earlier figure numbers |
| `workflows/` | Preparation, simulation, analysis and manuscript entry points |
| `ops/` | Cluster jobs, transfers and document utilities |
| `archive/` | Historical materials |

Raw trajectories and large simulation intermediates are excluded from Git. See [MD run data](docs/methods/md-run-data.md) for file formats and validation, and the [documentation index](docs/README.md) for methods and operations.
