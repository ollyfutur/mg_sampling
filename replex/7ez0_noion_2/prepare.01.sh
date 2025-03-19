#!/bin/bash

#forcefields="mMg_tip3p nMg_tip3p mMg_opc"
forcefields="ecc_spce"

for forcefield in $forcefields; do
   (
   cd $forcefield

   ### Generate PLUMED CVs include file
   # Bias scaling factors from optimization
   optimize_scalings rep??/bias.COLVAR.* \
      --bias_col cv \
      --scalings scalings.00.txt \
      -o scalings.01.txt \
      --start 25000 \
      --stride 5 > opt_scalings.01.log
   sh ../generate_plumed.sh scalings.01.txt > plumed.01.dat
   )
done
