#!/bin/bash

motifs="2r8s 3q3z 2r8s_noion 3q3z_noion"
#forcefields="mMg_tip3p"
forcefields="amber_tip3p"
from_rep=0
to_rep=3

for motif in $motifs; do
  mkdir $motif
  (
    cd $motif
    for forcefield in $forcefields; do
      mkdir $forcefield
      (
        cd $forcefield

        path="${PROJECT_DIR}/structures/motifs/${motif}/${forcefield}"

        ### Generate index file
        $GMX make_ndx -f ${path}/equil.tpr -o index.ndx << EOF
"RNA" & ! a H*
"RNA" & a O1P O2P
"RNA" & a O2 O4 O6
"SOL" & a OW
q
EOF
        sed -i "s/RNA_&_!H*/RNA_heavy/g" index.ndx
        sed -i "s/RNA_&_O1P_O2P/O_Phosphate/g" index.ndx
        sed -i "s/RNA_&_O2_O4_O6/O_Base/g" index.ndx
        sed -i "s/SOL_&_OW/O_Water/g" index.ndx

        ### Generate reference for PLUMED RMSD
        echo 'RNA_heavy' | gmx_mpi editconf -f ${path}/ions.gro -n index.ndx -o ref_heavy.pdb

        ### Generate PLUMED CVs include file
        cp ${PROJECT_DIR}/metad/op_${forcefield}/bias_params.wsv .
        # Bias scaling factors guess: sqrt
        seq $from_rep $to_rep | awk "{ print sqrt((NR-1)/($to_rep-$from_rep)) }" > scalings.txt
        sh ../../generate_plumed.sh > plumed.dat

        cp ${path}/processed.top .
        ### Generate each replica
        for i in $(seq -f "%02.0f" ${from_rep} ${to_rep}); do
          mkdir rep$i
          cp ${path}/equil.gro rep${i}/input.gro
        done
      )
    done
  )
done
