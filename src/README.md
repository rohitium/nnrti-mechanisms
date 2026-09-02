# `src/nnrti/`

The package. Importable as `nnrti` after `pip install -e .` or
`export PYTHONPATH=src`. Every entry point is
`python -m nnrti.<subpackage>.<module>`.

| Subpackage | Contents |
|---|---|
| `structure_prep/` | Structure and mutation preparation: run specifications, WT and mutant system building, manifest creation, mutation parsing and numbering |
| `md/` | MD runtime: task schema and CSV IO, the one-task execution entry point, and `openmm/` — staged heating and NPT production, minimisation, restraints, ligand parameterisation, MM/GBSA |
| `fep/` | Non-equilibrium alchemical free energy with pmx and GROMACS: hybrid topology preparation, switching, BAR/CGI/Jarzynski estimation, charge correction, panel assembly |
| `analysis/` | Analysis library: trajectory metrics, periodic-boundary handling, unit conversion, susceptibility loading, the genotype panel definition, result collection |
| `cli/` | The 16 scripts that produce a numbered figure or table. `paper/ARTIFACTS.md` maps each to its artifact |
| `utils/` | CIF parsing, mutation string handling, residue mapping, path helpers |

Directory constants are in `paths.py`.

## Trajectory time

Simulation time is derived from the run record rather than read from trajectory
metadata, which is not reliable for these files:

```python
production_ps = md_production_steps_completed * timestep_fs / 1000.0
time_ps = (frame_index + 1) * production_ps / n_frames
```

`analysis/md_timing.py` implements this; the profile workers in
`analysis/result_collector.py` use it.
