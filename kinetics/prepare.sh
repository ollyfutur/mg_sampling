#forcefields="amber_tip3p mMg_tip3p nMg_tip3p mMg_opc"
forcefields="ecc_spce"
otypes="op"
scalings="0.5 0.8 0.9 1.0 1.1 1.2 1.5"

# Needed to rescale coordination number to integer values
coef=1.039661

for ff in $forcefields; do
  for ox in $otypes; do
    dirname="${ox}_${ff}"
    mkdir -p $dirname
    cp $PROJECT_DIR/structures/diuridine/$ff/equil.gro $dirname/
    cp $PROJECT_DIR/structures/diuridine/$ff/processed.top $dirname/

    # Generate index file
    $GMX make_ndx -f $dirname/equil.gro -o $dirname/index.ndx << EOF
ri 2 & a O1P
ri 1 & a O4
"SOL" & a OW
q
EOF

    sed -i "s/r_2_&_O1P/O_Phosphate_1/g" $dirname/index.ndx &&
    sed -i "s/r_1_&_O4/O_Base_1/g" $dirname/index.ndx &&
    sed -i "s/SOL_&_OW/O_Water/g" $dirname/index.ndx &&

    # Parameters for the cosine envelope bias (from optimization)
    # Using parameters fitted for OP for all oxygen types
    cp $PROJECT_DIR/metad/op_${ff}/bias_params.wsv $dirname/
    read -r E sX sY <<< $(tail -n1 $dirname/bias_params.wsv)

    # Generate plumed file for each scaling value
    for scaling in $scalings; do
      mkdir $dirname/$scaling || continue
      cat >> $dirname/$scaling/plumed.dat << EOF
# Plumed file for testing bias kinetics
LOAD FILE=$PROJECT_DIR/ligandexchange/LigandExchange.so

mg: GROUP NDX_FILE=../index.ndx NDX_GROUP=MG
# First OP atom in topology
op1: GROUP NDX_FILE=../index.ndx NDX_GROUP=O_Phosphate_1
# First OB atom in topology
ob1: GROUP NDX_FILE=../index.ndx NDX_GROUP=O_Base_1
ow: GROUP NDX_FILE=../index.ndx NDX_GROUP=O_Water

d${ox}1: DISTANCE ATOMS=mg,${ox}1

cv: LIGANDEXCHANGE ...
   GROUPA=mg
   GROUPB=${ox}1
   GROUPC=ow
   SWITCH={COSINUS R_0=0.24 D_0=0.18}
   VAR=_nox,_now
   FUNC=-$E*step(Y)*step(0.5-Y)*step(X-0.00001)*(1-cos(Y*pi*4))*(1-cos(X*pi))/4;X=x+$sX*(1-cos(x*pi));Y=y-$sY*(1-cos(x*pi));x=$coef*_nox-$coef*_now+6;y=$coef*_nox+$coef*_now-6
   NLIST
   NL_CUTOFF_AB=1.0
   NL_CUTOFF_AC=1.0
   NL_STRIDE=500
   NESTED
   COMPONENTS
...

scaledcv: COMBINE ARG=cv.tot COEFFICIENTS=$scaling PERIODIC=NO

uwall: UPPER_WALLS ARG=d${ox}1 AT=0.8 KAPPA=50
bias: BIASVALUE ARG=scaledcv

#####################################################################
### Output ##########################################################
#####################################################################

PRINT ARG=d${ox}1,(cv\.n.*),(cv\.func.*),cv.tot,scaledcv,uwall.bias,bias.bias STRIDE=500 FILE=COLVAR
EOF
    done
  done
done
