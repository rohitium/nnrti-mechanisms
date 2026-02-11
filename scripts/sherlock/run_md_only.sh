#!/bin/bash
#
# Submit ONLY MD jobs on Sherlock using Snakemake
# Assumes prep was done locally and system XMLs already exist
#

set -e

echo "Submitting MD jobs only (skipping prep)..."
echo ""

# Target collect_and_analyze - it depends on all MD jobs completing
# Snakemake will submit only the missing MD jobs
snakemake collect_and_analyze \
  --profile workflow/profiles/sherlock \
  --rerun-incomplete \
  --until run_md

echo ""
echo "All MD jobs submitted. Monitor with:"
echo "  squeue -u \$USER"
echo "  tail -f .snakemake/slurm_logs/rule_run_md/*/*/*.log"
