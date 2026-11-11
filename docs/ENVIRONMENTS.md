# Environments

## Local work: nnrti-prep

Use **nnrti-prep** for structure and ligand preparation, trajectory analysis,
MM/GBSA, figures, tables and tests. It is the default environment for this project.

```bash
conda env create -f environment.yml          # first setup only
conda activate nnrti-prep
export PYTHONNOUSERSITE=1
python -m pip install --no-deps -e .         # register this source checkout
./workflows/05_manuscript_artifacts.sh
```

`environment.yml` specifies the local scientific software stack, including
OpenMM, OpenFF, AmberTools and trajectory libraries. `pyproject.toml` defines the
installable `nnrti` package and its Python dependencies. Keep both: creating the
environment and installing the project are different operations. `--no-deps`
keeps pip from replacing Conda-managed numerical packages. User-site packages
are disabled so `~/.local` cannot silently override the activated environment.
The environment recipe configures activation to disable user-site overrides;
the paper workflow also applies that isolation automatically.

The generated `paper/generated/reproduction.json` records the versions and hashes
actually used for each build. The environment specification is a supported minor
version recipe, not a complete cross-platform binary lock.

## Sherlock: site modules for simulation

Do not create a Conda environment on Sherlock for the normal workflow. Prepare
systems locally, transfer the serialized systems and matching topologies, and
run simulations in Slurm allocations using Sherlock modules. Transfer outputs
back for local analysis and publication builds.

The equilibrium MD submission scripts load (verified by importing both OpenMM
and the project simulation worker on Sherlock):

```bash
module load chemistry py-openmm/8.1.1_py312
export PYTHONNOUSERSITE=1
export PYTHONPATH="$GROUP_HOME/nnrti-mechanisms/code/src:${PYTHONPATH:-}"
python3 -c 'import openmm, nnrti; print(openmm.__version__)'
```

Preserve the module-provided `PYTHONPATH`: OpenMM is supplied through that path.
Submit from the scratch workspace, with `PROJECT_ROOT` explicitly set to it:

```bash
cd "$SCRATCH/nnrti-mechanisms-work"
export PROJECT_ROOT="$PWD"
# Holo: bash ops/slurm/cluster/submit_md_batched.sh
# Apo:  bash ops/slurm/cluster/submit_apo_md_groupsrc.sh
```

The workspace exposes the durable code through links while run outputs and logs
stay on scratch. Do not submit with Group Home as `PROJECT_ROOT`.

Prepared system XML/PDB assets let the simulation worker run without installing
local preparation dependencies on the cluster. Local OpenMM and the cluster
module have different versions; validate serialization and restart compatibility
before a new simulation campaign. Binary checkpoints are runtime-dependent.

For nonequilibrium FEP production, use the existing GROMACS module loader in a
separate shell or job:

```bash
source ops/slurm/cluster/load_gromacs_module.sh
```

It selects `chemistry gromacs/2023.1`, using `gmx` for preparation and `gmx_cuda`
for GPU production. Discover available versions with `ml spider` before changing
module selections. Do not combine the OpenMM and GROMACS CUDA stacks in one job.

## Optional local legacy tools: pmx

The legacy pmx hybrid-topology and raw FEP estimation tools use the separate
**pmx** environment (Python 3.11, NumPy < 2). This is a compatibility exception,
not a second environment for routine analysis. The paper build reads deposited
FEP estimates and integrated work values and does not require pmx.

`ops/slurm/fep/setup_pmx_env.sh` configures this local environment from
[deGrootLab/pmx](https://github.com/deGrootLab/pmx); the similarly named PyPI package
is unrelated. Use its pinned source revision and numerical versions and validate the
estimators before raw FEP reanalysis. Activating an old environment or successfully
importing pmx alone does not establish estimator compatibility.

The validated local compatibility stack is pmx `0+untagged.896.ge7b6a00`,
Python 3.11, NumPy 1.26.4, SciPy 1.13.1 and Matplotlib 3.8.4. It reproduces all
114 saved phase/replicate BAR estimates from integrated work within their reported
rounding. This check does not require the original per-switch trajectories.

PyMOL is an external optional tool for regenerating structure artwork; exact
submitted artwork is always available through the paper export. LibreOffice is
an optional external spreadsheet engine for native formula recalculation tests.

Historical cluster setup helpers and other installed environment names are not
part of the default workflow. Preserve their provenance before retiring them;
do not switch environments to resolve an error without identifying the missing
component first. Current storage paths and private inventories belong in the
durable project's `operations/` directory.
