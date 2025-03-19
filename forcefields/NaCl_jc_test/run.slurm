#!/bin/bash
#SBATCH --job-name=NaCljc
#SBATCH -N 1
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --ntasks-per-node=1
#SBATCH --gres=gpu:1
#SBATCH --account=EUHPC_R02_136
#SBATCH --partition=boost_usr_prod

module purge
module load gromacs/2023/2023.5-plumed-2.10b
export GMX_ENABLE_DIRECT_GPU_COMM=true

if [ ! -f run.tpr ]; then
  gmx_mpi grompp -f $PROJECT_DIR/forcefields/NaCl_jc_test/run.mdp -c ions.gro -o run.tpr
fi

srun gmx_mpi mdrun -deffnm run -cpi run.cpt -maxh 23.9
