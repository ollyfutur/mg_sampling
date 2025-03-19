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
      --scalings scalings.01.txt \
      -o scalings.02.txt \
      --start 15000 \
      --stride 5 > opt_scalings.02.log
   sh ../generate_plumed.sh scalings.02.txt > plumed.02.dat
   )
done
