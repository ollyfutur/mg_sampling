#!/bin/bash
#SBATCH --job-name=equil
#SBATCH -N 1
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --ntasks-per-node=4
#SBATCH --gres=gpu:4
#SBATCH --account=EUHPC_R02_136
#SBATCH --partition=boost_usr_prod

module purge
module load gromacs/2023/2023.5-plumed-2.10b
export GMX_ENABLE_DIRECT_GPU_COMM=true

for rep in rep??; do
  if [ ! -f $rep/equil.tpr ]; then
    gmx_mpi grompp -f $PROJECT_DIR/common/equil_posre.mdp -p topol.top -c $rep/min.gro -r $rep/min.gro -o $rep/equil.tpr
  fi
  if [ ! -f $rep/equil.gro ]; then
    (
      cd $rep
      srun gmx_mpi mdrun -deffnm equil
    )
  fi
done

