systems="7ez0 7ez0_noion 7ez0_2 7ez0_noion_2"
forcefields="mMg_opc mMg_tip3p nMg_tip3p ecc_spce"

for system in $systems; do
  for forcefield in $forcefields; do
    dir="$system/$forcefield"
    plumedfile=$(echo $dir/plumed.??.dat | awk '{ print $1 }')
    sed -e '/NL/d' -e 's/NESTED/COMPONENTS/' -e '/rmsd/d' -e '/scaledcv/d' -e '/bias/d' $plumedfile > $dir/plumed.analysis.dat
    echo "PRINT ARG=(cv\.nb.*),(cv\.nc.*) STRIDE=5000 FILE=cv.COLVAR" >> $dir/plumed.analysis.dat
  done
done
