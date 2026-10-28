#!/bin/bash

system="1naj"
forcefield="mMg_tip3p"
from_rep=0
to_rep=2
n_rep=$((to_rep-from_rep+1))
# Intensity of bias
lambda_index=7

mkdir -p $system
path=$PROJECT_DIR/structures/$system/$forcefield
mkdir -p $system/${forcefield}
(
   cd $system/${forcefield}

   ### Generate index file
   gmx_mpi make_ndx -f $path/rep00/equil.tpr -o index.ndx << EOF
"DNA" & a P C1' C2
"DNA" & a O1P O2P
"DNA" & a O2 O4 O6
"SOL" & a OW
q
EOF

   sed -i "s/DNA_&_P_C1'_C2/DNA_subset/g" index.ndx
   sed -i "s/DNA_&_O1P_O2P/O_Phosphate/g" index.ndx
   sed -i "s/DNA_&_O2_O4_O6/O_Base/g" index.ndx
   sed -i "s/SOL_&_OW/O_Water/g" index.ndx

   # Generate plumed file
   sh ../../generate_plumed.sh \
      $PROJECT_DIR/metad/op_${forcefield}/bias_params.wsv \
      $PROJECT_DIR/kinetics/op_${forcefield}/opt_scalings_ramp.txt \
      $PROJECT_DIR/kinetics/op_${forcefield}/opt_alpha_ramp.txt \
      $lambda_index > plumed.dat

   cp $path/processed.top .
   ### Generate each replicate
   for i in $(seq -f "%02.0f" ${from_rep} ${to_rep}); do
      mkdir -p rep$i
      cp $path/rep$i/equil.gro rep$i/input.gro
   done
)
