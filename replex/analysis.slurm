#!/bin/bash
#SBATCH --job-name=anal
#SBATCH -N 1
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=2
#SBATCH --ntasks-per-node=12
#SBATCH --account=IscrB_MGReplex
#SBATCH --partition=boost_usr_prod

module purge
module load gromacs/2023/2023.5-plumed-2.10b
export PLUMED_NUM_THREADS=2

sh $PROJECT_DIR/replex/generate_plumed_analysis.sh bias_params.wsv > plumed_analysis.dat
for rep in rep??; do
   cd $rep
   srun -n 1 -c $PLUMED_NUM_THREADS --exclusive plumed driver --plumed ../plumed_analysis.dat --ixtc run.xtc --trajectory-stride 0 --timestep 0.002 &
   cd ..
done

wait

# Same suffix as bias.COLVAR.<i>: index of the replica within the replicate
i=0
for rep in rep??; do
   mv $rep/cv.COLVAR $rep/cv.COLVAR.${i}
   i=$((i+1))
done
