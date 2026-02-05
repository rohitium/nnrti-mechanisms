#!/bin/bash
#SBATCH --job-name=v106a_rep1
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --time=4:00:00
#SBATCH --mem=16G
#SBATCH --output=logs/v106a_rep1_%j.out
#SBATCH --error=logs/v106a_rep1_%j.err

ml chemistry py-openmm/8.1.1_py312
export OPENMM_PLATFORM=CUDA
mkdir -p logs

# V106A replicate 1 uses task IDs 6 (complex) and 7 (solvent)
python3 -m src.cluster.fep_worker \
  --manifest results/fep_manifest.csv \
  --task-id 6 \
  --equil-steps 10000 \
  --prod-steps 25000 \
  --sample-interval 200

python3 -m src.cluster.fep_worker \
  --manifest results/fep_manifest.csv \
  --task-id 7 \
  --equil-steps 10000 \
  --prod-steps 25000 \
  --sample-interval 200
