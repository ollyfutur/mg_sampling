#!/bin/bash
#SBATCH --job-name=anal
#SBATCH -N 1
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=4
#SBATCH --ntasks-per-node=8
#SBATCH --account=Sis25_bussi
#SBATCH --partition=boost_usr_prod
#SBATCH --qos=boost_qos_dbg

module purge
module load gromacs/2023/2023.5-plumed-2.10b

cp $PROJECT_DIR/replex/7ez0/mMg_tip3p/replicate0/plumed_analysis.dat .

export PLUMED_NUM_THREADS=4
for rep in rep??; do
   cd $rep
   srun -n 1 --exclusive plumed driver --plumed ../plumed_analysis.dat --ixtc run.xtc --trajectory-stride 0 --timestep 0.002 &
   cd ..
done

wait
