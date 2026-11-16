#!/bin/bash
# One replicate of a system per array task, each on 1 GPU (array index = replicate number).
# Usage: sbatch --array=0-7 density_maps.slurm plain_md/7ez0/mMg_tip3p
# All systems: submit_density.sh. Existing outputs are skipped by the script.
#SBATCH --job-name=density
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --mem=120G
#SBATCH --time=02:00:00
#SBATCH --account=IscrB_MGReplex
#SBATCH --partition=boost_usr_prod
#SBATCH --output=density_%A_%a.out

set -euo pipefail

# Adapt to the Python environment on the cluster
source "$PROJECT_DIR/.venv/bin/activate"

# One thread per worker process (avoid oversubscription)
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1

cd "$PROJECT_DIR/density"
python density_maps.py "$1" "${SLURM_ARRAY_TASK_ID}" --n-workers "${SLURM_CPUS_PER_TASK}"
