#!/bin/bash
#
# Submit ONLY MD jobs on Sherlock using Snakemake
# Assumes prep was done locally and system XMLs already exist
#

set -e

echo "Submitting MD jobs only (skipping prep)..."
echo ""

# Just target the run_md rule - Snakemake will find all incomplete MD jobs
snakemake run_md \
  --profile workflow/profiles/sherlock \
  --rerun-incomplete

echo ""
echo "All MD jobs submitted. Monitor with:"
echo "  squeue -u \$USER"
echo "  tail -f .snakemake/slurm_logs/rule_run_md/*/*/*.log"
