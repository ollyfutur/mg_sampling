#!/bin/bash
#SBATCH --job-name=anal
#SBATCH -N 1
#SBATCH --time=00:30:00
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=2
#SBATCH --ntasks-per-node=12
#SBATCH --account=Sis25_bussi
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

for rep in rep??; do
   i=$((10#${rep:(-2)}))
   mv $rep/cv.COLVAR $rep/cv.COLVAR.${i}
done
