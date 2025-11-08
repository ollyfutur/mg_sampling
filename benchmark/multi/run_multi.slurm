#!/bin/bash
#SBATCH --job-name=plainmd
#SBATCH -N 2
#SBATCH --time=24:00:00
#SBATCH --time=00:30:00
#SBATCH --gres=gpu:4
#SBATCH --cpus-per-gpu=8
#SBATCH --account=Sis25_bussi
#SBATCH --partition=boost_usr_prod
#SBATCH --qos=boost_qos_dbg

module purge
module load gromacs/2023/2023.5-plumed-2.10b
export GMX_ENABLE_DIRECT_GPU_COMM=true

# srun -n 8 gmx_mpi mdrun -deffnm run -multidir rep? -plumed ../plumed.dat -nsteps 20000 -resetstep 15000
# cp rep0/run.log n8.log
# srun -n 8 gmx_mpi mdrun -deffnm run -multidir rep? -plumed ../plumed.dat -npme 1 -ntomp_pme 8 -ntomp 8 -nsteps 20000 -resetstep 15000
# cp rep0/run.log n8_npme1_cpme8.log
srun -n 16 gmx_mpi mdrun -deffnm run -multidir rep? -plumed ../plumed.dat -npme 1 -ntomp_pme 1 -ntomp 9 -nsteps 20000 -resetstep 15000
cp rep0/run.log n8_npme1_cpme1.log
