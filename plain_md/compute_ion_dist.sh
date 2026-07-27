ions="MG NA"
otypes="O_Phosphate"

replicas="rep??"
for rep in $replicas; do
  for otype in $otypes; do
    for ion in $ions; do
      $GMX pairdist -f $rep/run.xtc -s $rep/run.tpr -n $rep/../index.ndx -cutoff 2 -refgrouping none -selgrouping all -ref $otype -sel $ion -o $rep/dist_${otype}_${ion}.xvg &
      $GMX pairdist -f $rep/run.xtc -s $rep/run.tpr -n $rep/../index.ndx -cutoff 2 -refgrouping none -selgrouping all -ref $ion -sel $otype -o $rep/dist_${ion}_${otype}.xvg &
    done
  done
done
