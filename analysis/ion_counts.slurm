#!/bin/bash
# One replicate of a system per array task (array index = replicate number).
# Usage: sbatch --array=0-7 ion_counts.slurm plain_md/7ez0/mMg_tip3p
# All systems: submit_ion_counts.sh. Existing outputs are skipped by the script.
#SBATCH --job-name=ioncount
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --account=IscrB_MGReplex
#SBATCH --partition=boost_usr_prod
#SBATCH --output=ioncount_%A_%a.out

set -euo pipefail

# Adapt to the Python environment on the cluster
source "$PROJECT_DIR/.venv/bin/activate"

cd "$PROJECT_DIR/analysis"
python ion_counts.py "$1" "${SLURM_ARRAY_TASK_ID}"
