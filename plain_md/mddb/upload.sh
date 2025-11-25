for rep in rep00 rep01 rep02 rep03; do
  sh curl.sh 7ez0_noion/nMg_tip3p/rep00/solute.tpr ${rep}.tpr
  sh curl.sh 7ez0_noion/nMg_tip3p/rep00/wrapped.xtc ${rep}.xtc
  sh curl.sh 7ez0_noion/nMg_tip3p/rep00/rmsd.COLVAR.? ${rep}.COLVAR
done
