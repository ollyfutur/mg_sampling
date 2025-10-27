#!/bin/bash
#SBATCH --job-name=metad_static
#SBATCH -N 1
#SBATCH --time=00:30:00
##SBATCH --cpus-per-task=4
#SBATCH --cpus-per-task=32
##SBATCH --ntasks-per-node=8
#SBATCH --ntasks-per-node=1
#SBATCH --account=Sis25_bussi
#SBATCH --partition=boost_usr_prod
#SBATCH --qos=boost_qos_dbg

module purge
module load gromacs/2023/2023.5-plumed-2.10b
export GMX_ENABLE_DIRECT_GPU_COMM=true

#export PLUMED_NUM_THREADS=4
export PLUMED_NUM_THREADS=32

#dirnames="ob_mMg_opc ob_nMg_tip3p op_mMg_opc op_nMg_tip3p ob_amber_tip3p ob_mMg_tip3p op_amber_tip3p op_mMg_tip3p"
#dirnames="ob_mMg_tip3p"
dirnames="op_mMg_tip3p"
for dirname in $dirnames; do
  echo $dirname
  (
  cd $dirname
  srun --ntasks=1 --cpus-per-task=$PLUMED_NUM_THREADS plumed driver --plumed plumed_static.dat --mf_xtc run.xtc --timestep 0.002 --trajectory-stride 500 > plumed_static.log
  ) &
done
wait
