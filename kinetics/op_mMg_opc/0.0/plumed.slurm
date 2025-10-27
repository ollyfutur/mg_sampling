#!/bin/bash
#SBATCH --job-name=kinetics
#SBATCH -N 1
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=8
#SBATCH --ntasks-per-node=1
#SBATCH --gres=gpu:1
#SBATCH --account=EUHPC_R02_136
#SBATCH --partition=boost_usr_prod
#SBATCH --qos=boost_qos_dbg

module purge
module load gromacs/2023/2023.5-plumed-2.10b
export GMX_ENABLE_DIRECT_GPU_COMM=true

plumed driver --ixtc run.xtc --plumed plumed_analysis.dat --trajectory-stride 0 --timestep 0.002
