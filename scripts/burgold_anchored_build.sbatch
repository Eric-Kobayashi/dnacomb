#!/usr/bin/env bash
#SBATCH --job-name=dnacomb-anchored-build
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=logs/anchored-build-%j.log
set -euo pipefail
analysis="$HOME/lustre/dnacomb/analysis/burgold"
project="$HOME/lustre/dnacomb/software/dnacomb"
pixi="$HOME/lustre/dnacomb/software/pixi/bin/pixi"
trap 'status=$?; printf "%s\n" "$status" > "$analysis/anchored_build.exitcode"; exit "$status"' EXIT
export CARGO_BUILD_JOBS="$SLURM_CPUS_PER_TASK"
cd "$project/source"
"$pixi" run --manifest-path "$project/pixi.toml" --locked cargo fmt --check
"$pixi" run --manifest-path "$project/pixi.toml" --locked cargo test --locked
"$pixi" run --manifest-path "$project/pixi.toml" --locked build
