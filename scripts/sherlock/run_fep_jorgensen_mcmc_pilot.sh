#!/bin/bash
#
# MCMC is NOT supported on Sherlock with the py-openmm module alone.
#
# Perses HybridRepexSampler needs perses + openmmtools + ambertools (nnrti-prep stack).
# Sherlock discourages conda; our validated GPU path is fixed-λ worker.py instead.
#
# Run MCMC locally (Mac with nnrti-prep + GPU):
#   PYTHONNOUSERSITE=1 PYTHONPATH=. python -m scripts.fep_jorgensen.check_mcmc_env
#   PYTHONNOUSERSITE=1 PYTHONPATH=. python -m scripts.fep_jorgensen.mcmc_sample \
#     --mutation V106A --n-cycles 20 --steps-per-cycle 50
#
echo "ERROR: MCMC (mcmc_sample.py) is not supported on Sherlock with module load py-openmm." >&2
echo "Use worker.py fixed-lambda windows on Sherlock, or run mcmc_sample.py locally." >&2
exit 1
