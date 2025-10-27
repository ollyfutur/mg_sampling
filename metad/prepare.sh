#forcefields="amber_tip3p mMg_tip3p nMg_tip3p mMg_opc"
#forcefields="ecc_spce"
#otypes="op ob"
forcefields="mMg_tip3p"
otypes="op"

for ff in $forcefields; do
  for ox in $otypes; do
    dirname="${ox}_${ff}"
    mkdir $dirname || continue

    cp $PROJECT_DIR/structures/diuridine/$ff/equil.gro $dirname/
    cp $PROJECT_DIR/structures/diuridine/$ff/processed.top $dirname/

    # Generate index file
    $GMX make_ndx -f $dirname/equil.gro -o $dirname/index.ndx << EOF
"RNA" & a O1P O2P
"RNA" & a O2 O4 O6
"SOL" & a OW
q
EOF

    sed -i "s/RNA_&_O1P_O2P/O_Phosphate/g" $dirname/index.ndx &&
    sed -i "s/RNA_&_O2_O4_O6/O_Base/g" $dirname/index.ndx &&
    sed -i "s/SOL_&_OW/O_Water/g" $dirname/index.ndx &&

    # Generate plumed file
    cat >> $dirname/plumed.dat << EOF
# Plumed file for Metadynamics

mg: GROUP NDX_FILE=index.ndx NDX_GROUP=MG
op: GROUP NDX_FILE=index.ndx NDX_GROUP=O_Phosphate
ob: GROUP NDX_FILE=index.ndx NDX_GROUP=O_Base
ow: GROUP NDX_FILE=index.ndx NDX_GROUP=O_Water

_n${ox}: COORDINATION ...
   GROUPA=mg
   GROUPB=${ox}
   SWITCH={CUSTOM FUNC=0.5*(cos(pi*x/0.24)+1) R_0=1 D_0=0.18 D_MAX=0.42}
...

_now: COORDINATION ...
   GROUPA=mg
   GROUPB=ow
   SWITCH={CUSTOM FUNC=0.5*(cos(pi*x/0.24)+1) R_0=1 D_0=0.18 D_MAX=0.42}
...

# SWITCH={COSINUS R_0=0.24 D_0=0.18} <- wrong derivatives: https://github.com/plumed/plumed2/issues/1317
_n${ox}_old: COORDINATION ...
   GROUPA=mg
   GROUPB=${ox}
   SWITCH={COSINUS R_0=0.24 D_0=0.18}
...

_now_old: COORDINATION ...
   GROUPA=mg
   GROUPB=ow
   SWITCH={COSINUS R_0=0.24 D_0=0.18}
...


n${ox}: COMBINE ARG=_n${ox} COEFFICIENTS=1.039661 PERIODIC=NO
now: COMBINE ARG=_now COEFFICIENTS=1.039661 PERIODIC=NO

lwall: LOWER_WALLS ARG=now AT=4 KAPPA=500
uwall: UPPER_WALLS ARG=now AT=7 KAPPA=500

metad: METAD ...
   ARG=n${ox},now
   SIGMA=0.05,0.05 HEIGHT=1.0 PACE=500 TEMP=300 BIASFACTOR=15 
...

PRINT ARG=n${ox},now,lwall.bias,uwall.bias,metad.bias,_n${ox},_now,_n${ox}_old,_now_old STRIDE=500 FILE=COLVAR
EOF

  done
done
