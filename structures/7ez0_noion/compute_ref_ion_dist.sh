otypes="O_Phosphate O_Base"
ion="MG"

$GMX make_ndx -f rep00/ions1.tpr -o ref_index.ndx << EOF
"RNA" & a O1P O2P
"RNA" & a O2 O4 O6
q
EOF
sed -i "s/RNA_&_O1P_O2P/O_Phosphate/g" ref_index.ndx
sed -i "s/RNA_&_O2_O4_O6/O_Base/g" ref_index.ndx

for otype in $otypes; do
    $GMX pairdist -f box_solv.gro -s rep00/ions1.tpr -n ref_index.ndx -cutoff 2 -refgrouping none -selgrouping all -ref $otype -sel $ion -o ref_dist_${otype}_${ion}.xvg
done
