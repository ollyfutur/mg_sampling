#!/bin/bash
#SBATCH --job-name=min
#SBATCH -N 1
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=16
#SBATCH --account=EUHPC_R02_136
#SBATCH --partition=boost_usr_prod

module purge
module load gromacs/2023/2023.5-plumed-2.10b
export GMX_ENABLE_DIRECT_GPU_COMM=true

for rep in rep??; do
  if [ ! -f $rep/min.tpr ]; then
    if [ ! -f processed.top ]; then
      gmx_mpi grompp -f $PROJECT_DIR/common/min.mdp -p topol.top -pp processed.top -c $rep/ions2.gro -o $rep/min.tpr
    else
      gmx_mpi grompp -f $PROJECT_DIR/common/min.mdp -p topol.top -c $rep/ions2.gro -o $rep/min.tpr
    fi
  fi
done

srun -n 16 gmx_mpi mdrun -multidir rep?? -deffnm min
