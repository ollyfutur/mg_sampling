#!/bin/bash
#SBATCH --job-name=fwddens
#SBATCH -N 1
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --ntasks-per-node=1
#SBATCH --gres=gpu:1
#SBATCH --account=Sis25_bussi
#SBATCH --partition=boost_usr_prod
#SBATCH --time=00:30:00
#SBATCH --qos=boost_qos_dbg


#srun python forward_density_plain_md.py
srun python forward_density_replex.py
