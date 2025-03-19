forcefields="amber_tip3p mMg_tip3p nMg_tip3p mMg_opc ecc_spce"
otypes="op ob"

for ff in $forcefields; do
  for ox in $otypes; do
    dirname="${ox}_${ff}"
    sed -e '/METAD/a \   RESTART=YES' \
        -e 's/HEIGHT=1.0/HEIGHT=0.0/g' \
        -e 's/PACE=500/PACE=1000000000/g' \
        -e 's/COLVAR/static.COLVAR/g' \
        $dirname/plumed.dat > $dirname/plumed_static.dat
  done
done
