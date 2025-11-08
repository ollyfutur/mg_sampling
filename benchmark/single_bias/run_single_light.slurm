#!/bin/bash
#SBATCH --job-name=plainmd
#SBATCH -N 1
#SBATCH --time=24:00:00
#SBATCH --time=00:30:00
#SBATCH --gres=gpu:4
#SBATCH --cpus-per-gpu=8
#SBATCH --account=Sis25_bussi
#SBATCH --partition=boost_usr_prod
#SBATCH --qos=boost_qos_dbg

module purge
# module load plumed/2.10/2.10bnl-patch
module load plumed/2.10/2.10b
module load gromacs/2023/2023.5-plumed-2.10b
export GMX_ENABLE_DIRECT_GPU_COMM=true

# srun -n 4 -G 4 gmx_mpi mdrun -s run.tpr -plumed ../plumed_bias.dat -nsteps 25000 -resetstep 20000
# mv md.log n4.log
srun -n 1 -G 1 gmx_mpi mdrun -s run.tpr -plumed ../plumed_bias.dat -nsteps 25000 -resetstep 20000
mv md.log n1.log
# srun -n 4 -G 4 gmx_mpi mdrun -s run.tpr -plumed ../plumed_bias.dat -npme 1 -ntomp_pme 2 -ntomp 10 -nsteps 25000 -resetstep 20000
# mv md.log n4_pme1.log

