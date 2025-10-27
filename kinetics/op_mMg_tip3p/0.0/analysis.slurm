#!/bin/bash
#SBATCH --job-name=anal
#SBATCH -N 1
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --account=Sis25_bussi
#SBATCH --partition=boost_usr_prod
#SBATCH --qos=boost_qos_dbg

module purge
module load gromacs/2023/2023.5-plumed-2.10b

srun -n 1 --exclusive plumed driver --plumed plumed_analysis.dat --ixtc run.xtc --trajectory-stride 0 --timestep 0.002
