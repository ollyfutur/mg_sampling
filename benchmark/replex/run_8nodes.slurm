#!/bin/bash
#SBATCH --job-name=replex
#SBATCH -N 8
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --ntasks-per-node=4
#SBATCH --gres=gpu:4
#SBATCH --account=Sis25_bussi
#SBATCH --partition=boost_usr_prod

module purge
module load gromacs/2023/2023.5-plumed-2.10b
export GMX_ENABLE_DIRECT_GPU_COMM=true

srun gmx_mpi mdrun -multidir rep?? -deffnm run -cpi run.cpt -plumed ../plumed.dat -nsteps 20000 -resetstep 15000
cp rep00/run.log 8nodes.log
