#!/bin/bash

#forcefields="mMg_tip3p nMg_tip3p"
forcefields="ecc_spce"
#forcefields="mMg_opc"
from_rep=0
to_rep=11
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
"RNA" & a O2' O4' O3' O5'
"RNA" & a N*
"SOL" & a OW
q
EOF

      sed -i "s/RNA_&_P_C1'_C2/RNA_subset/g" index.ndx
      sed -i "s/RNA_&_O1P_O2P/O_Phosphate/g" index.ndx
      sed -i "s/RNA_&_O2_O4_O6/O_Base/g" index.ndx
      sed -i "s/RNA_&_O2'_O4'_O3'_O5'/O_Ribose/g" index.ndx
      sed -i "s/RNA_&_N\*/N_Base/g" index.ndx
      sed -i "s/SOL_&_OW/O_Water/g" index.ndx

      ### Generate reference for PLUMED RMSD
      echo 'RNA_subset' | gmx_mpi editconf -f $path/rep00/ions2.gro -n index.ndx -o ref_subset.pdb

      ### Generate PLUMED CVs include file
      cp $PROJECT_DIR/metad/op_${forcefield}/bias_params.wsv .
      # Bias scaling factors guess: sqrt
      seq $from_rep $to_rep | awk "{ print sqrt((NR-1)/($to_rep-$from_rep)) }" > scalings.00.txt
      sh ../generate_plumed.sh scalings.00.txt > plumed.00.dat

      cp $path/processed.top .
      ### Generate each replica
      for i in $(seq -f "%02.0f" ${from_rep} ${to_rep}); do
         mkdir rep$i
         cp $path/rep$i/equil.gro rep$i/input.gro
      done
   )
done
