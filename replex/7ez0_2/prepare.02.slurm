#!/bin/bash
#SBATCH --job-name=prep
#SBATCH -N 1
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --ntasks-per-node=1
#SBATCH --gres=gpu:1
#SBATCH --account=EUHPC_R02_136
#SBATCH --partition=boost_usr_prod

#forcefields="mMg_tip3p nMg_tip3p mMg_opc"
forcefields="ecc_spce"

for forcefield in $forcefields; do
   (
   cd $forcefield

   ### Generate PLUMED CVs include file
   # Bias scaling factors from optimization
   optimize_scalings rep??/bias.COLVAR.* \
      --bias_col cv \
      --scalings scalings.01.txt \
      -o scalings.02.txt \
      --start 100000 \
      --stride 5 > opt_scalings.02.log
   sh ../generate_plumed.sh scalings.02.txt > plumed.02.dat
   )
done
