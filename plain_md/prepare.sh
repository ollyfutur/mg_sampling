#!/bin/bash

#systems="7ez0 7ez0_noion"
systems="7ez0_noion"
forcefields="nMg_tip3p"
from_rep=0
to_rep=3
n_rep=$((to_rep-from_rep+1))

for system in $systems; do
   mkdir -p $system
   for forcefield in $forcefields; do
      path=$PROJECT_DIR/structures/$system/$forcefield
      mkdir -p $system/${forcefield}
      (
         cd $system/${forcefield}

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

         sh ../../generate_plumed.sh > plumed.dat

         cp $path/processed.top .
         ### Generate each replica
         for i in $(seq -f "%02.0f" ${from_rep} ${to_rep}); do
            mkdir -p rep$i
            cp $path/rep$i/equil.gro rep$i/input.gro
         done
      )
   done
done
