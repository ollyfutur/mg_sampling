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

for rep in rep??; do
    gmx_mpi grompp -f $PROJECT_DIR/common/run_light.mdp -p processed.top -c $rep/input.gro -o $rep/run.tpr
    gmx_mpi convert-tpr -s $rep/run.tpr -until 200000 -o $rep/run.tpr
done

srun gmx_mpi mdrun -multidir rep?? -deffnm run -cpi run.cpt -plumed ../plumed.dat -nb gpu -pme gpu -npme 1 -ntomp 10 -ntomp_pme 2 -maxh 23.9
