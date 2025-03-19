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
      --scalings scalings.02.txt \
      -o scalings.03.txt \
      --start 50000 \
      --stride 5 > opt_scalings.03.log
   sh ../generate_plumed.sh scalings.03.txt > plumed.03.dat
   )
done
