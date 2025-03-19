systems="2r8s 2r8s_noion 3q3z 3q3z_noion"
forcefields="mMg_tip3p"
ions="MG CL"
otypes="O_Phosphate O_Base"

for system in $systems; do
  for ff in $forcefields; do
    replicas="$system/$ff/rep??"
    for rep in $replicas; do
      for otype in $otypes; do
        for ion in $ions; do
          $GMX pairdist -f $rep/run.xtc -s $rep/run.tpr -n $rep/../index.ndx -cutoff 2 -refgrouping none -selgrouping all -ref $otype -sel $ion -o $rep/dist_${otype}_${ion}.xvg > /dev/null &
        done
      done
    done
  done
done
wait

