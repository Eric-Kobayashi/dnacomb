#!/usr/bin/env bash
#SBATCH --job-name=burgold-anchored
#SBATCH --partition=cpu
#SBATCH --array=0-28%3
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=00:30:00
#SBATCH --output=logs/anchored-%A_%a.log
set -euo pipefail
cd "$HOME/lustre/dnacomb/analysis/burgold"
python3 run_anchored.py --index "$SLURM_ARRAY_TASK_ID"
