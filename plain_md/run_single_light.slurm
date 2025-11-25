#!/bin/bash
#SBATCH --job-name=plainmd
#SBATCH -N 1
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --ntasks-per-node=1
#SBATCH --gres=gpu:1
#SBATCH --account=Sis25_bussi
#SBATCH --partition=boost_usr_prod

module purge
module load gromacs/2023/2023.5-plumed-2.10b
export GMX_ENABLE_DIRECT_GPU_COMM=true

# Use lighter output xtc frequency (100ps)

if [ ! -f run.tpr ]; then
  gmx_mpi grompp -f $PROJECT_DIR/common/run_light.mdp -p ../processed.top -c input.gro -o run.tpr
fi
gmx_mpi convert-tpr -s run.tpr -until 600000 -o run.tpr

srun gmx_mpi mdrun -deffnm run -cpi run.cpt -plumed ../plumed.dat -maxh 23.9
