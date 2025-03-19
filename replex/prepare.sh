#!/bin/bash

#forcefields="mMg_tip3p nMg_tip3p mMg_opc"
forcefields="mMg_tip3p_noion nMg_tip3p_noion mMg_opc_noion"
from_rep=0
to_rep=7
n_rep=$((to_rep-from_rep+1))

for forcefield in $forcefields; do
   path=$PROJECT_DIR/structures/7ez0/$forcefield
   mkdir $forcefield
   (
      cd $forcefield

      ### Generate index file
      gmx_mpi make_ndx -f $path/rep00/equil.tpr -o index.ndx << EOF
"RNA" & a P C1' C2
"RNA" & a O1P O2P
"RNA" & a O2 O4 O6
"SOL" & a OW
q
EOF

      sed -i "s/RNA_&_P_C1'_C2/RNA_subset/g" index.ndx
      sed -i "s/RNA_&_O1P_O2P/O_Phosphate/g" index.ndx
      sed -i "s/RNA_&_O2_O4_O6/O_Base/g" index.ndx
      sed -i "s/SOL_&_OW/O_Water/g" index.ndx

      ### Generate reference for PLUMED RMSD
      echo 'RNA_subset' | gmx_mpi editconf -f $path/rep00/ions2.gro -n index.ndx -o ref_subset.pdb

      ### Generate PLUMED CVs include file
      cp $PROJECT_DIR/metad/op_${forcefield%_noion}/bias_params.wsv .
      # Bias scaling factors from optim
      cp $PROJECT_DIR/replex/7ez0/optimize/${n_rep}_replicas/${forcefield%_noion}/opt_scalings.txt scalings.txt
      sh ../generate_plumed.sh > plumed.dat

      cp $path/processed.top .
      ### Generate each replica
      for i in $(seq -f "%02.0f" ${from_rep} ${to_rep}); do
         mkdir rep$i
         cp $path/rep$i/equil.gro rep$i/input.gro
      done
   )
done
