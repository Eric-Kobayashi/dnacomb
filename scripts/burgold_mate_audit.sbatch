#!/usr/bin/env bash
#SBATCH --job-name=burgold-mate-audit
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:10:00
#SBATCH --output=logs/anchored-audit-%j.log
set -euo pipefail
cd "$HOME/lustre/dnacomb/analysis/burgold"
trap 'status=$?; printf "%s\n" "$status" > results_anchored/coordinate_audit.exitcode; exit "$status"' EXIT
python3 audit_burgold_mate_sequences.py
