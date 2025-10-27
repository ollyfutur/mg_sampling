#!/bin/bash
#SBATCH --job-name=kinetics
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

if [ ! -f run.tpr ]; then
  gmx_mpi grompp -f $PROJECT_DIR/kinetics/run_bench.mdp -p ../processed.top -c ../equil.gro -o run.tpr -maxwarn 1
fi

srun gmx_mpi mdrun -deffnm run -cpi run.cpt -plumed plumed.dat -maxh 0.5
