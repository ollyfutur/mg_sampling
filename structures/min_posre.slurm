#!/bin/bash
#SBATCH --job-name=min
#SBATCH -N 1
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=8
#SBATCH --ntasks-per-node=1
#SBATCH --gres=gpu:1
#SBATCH --account=Sis25_bussi
#SBATCH --partition=boost_usr_prod
#SBATCH --qos=boost_qos_dbg

module purge
module load gromacs/2023/2023.5-plumed-2.10b
export GMX_ENABLE_DIRECT_GPU_COMM=true

for rep in rep??; do
  if [ ! -f $rep/min.tpr ]; then
    gmx_mpi grompp -f $PROJECT_DIR/common/min_posre.mdp -p topol.top -c $rep/ions2.gro -r $rep/ions2.gro -o $rep/min.tpr
  fi
  if [ ! -f $rep/min.gro ]; then
    (
      cd $rep
      gmx_mpi mdrun -deffnm min
    )
  fi
done

