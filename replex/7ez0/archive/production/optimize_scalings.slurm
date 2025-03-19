#!/bin/bash
#SBATCH --job-name=opt
#SBATCH -N 1
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=16
#SBATCH --ntasks-per-node=1
#SBATCH --gres=gpu:0
#SBATCH --account=EUHPC_R02_136
#SBATCH --partition=boost_usr_prod

source $PROJECT_DIR/.venv/bin/activate

python ../optimize_scalings.py > optimize_scalings.log
