#forcefields="amber_tip3p mMg_tip3p nMg_tip3p mMg_opc"
forcefields="mMg_tip3p"
scalings="1.0"
strides="100 250 500 1000"
cutoffs="1.1 1.0 0.9 0.8 0.7"

# Needed to rescale coordination number to integer values
coef=1.039661

for ff in $forcefields; do
  dirname="benchmark_nlist_op_${ff}"
  mkdir -p $dirname
  cp $PROJECT_DIR/structures/diuridine/$ff/equil.gro $dirname/
  cp $PROJECT_DIR/structures/diuridine/$ff/processed.top $dirname/

  # Generate index file
  $GMX make_ndx -f $dirname/equil.gro -o $dirname/index.ndx <<EOF
ri 2 & a O1P O2P
"SOL" & a OW
q
EOF

  sed -i "s/r_2_&_O1P_O2P/O_Phosphate/g" $dirname/index.ndx &&
  sed -i "s/SOL_&_OW/O_Water/g" $dirname/index.ndx &&

  # Parameters for the cosine envelope bias (from optimization)
  # Using parameters fitted for OP for all oxygen types
  cp $PROJECT_DIR/metad/op_${ff}/bias_params.wsv $dirname/
  read -r E sX sY <<<$(tail -n1 $dirname/bias_params.wsv)

  # Generate plumed file for each scaling value
  for scaling in $scalings; do
    mkdir $dirname/$scaling || continue
    cat >>$dirname/$scaling/plumed.dat <<EOF
# Plumed file for testing bias kinetics
LOAD FILE=$PROJECT_DIR/ligandexchange/LigandExchange.so

mg: GROUP NDX_FILE=../index.ndx NDX_GROUP=MG
# First OP atom in topology
op: GROUP NDX_FILE=../index.ndx NDX_GROUP=O_Phosphate
ow: GROUP NDX_FILE=../index.ndx NDX_GROUP=O_Water

op_com: COM ATOMS=op
dop: DISTANCE ATOMS=mg,op_com

# Reference
cv: LIGANDEXCHANGE ...
   GROUPA=mg
   GROUPB=op
   GROUPC=ow
   SWITCH={CUSTOM FUNC=0.5*(cos(pi*x/0.24)+1) R_0=1 D_0=0.18 D_MAX=0.42}
   VAR=_nop,_now
   FUNC=-$E*step(Y)*step(0.5-Y)*step(X-0.00001)*(1-cos(Y*pi*4))*(1-cos(X*pi))/4;X=x+$sX*(1-cos(x*pi));Y=y-$sY*(1-cos(x*pi));x=$coef*_nop-$coef*_now+6;y=$coef*_nop+$coef*_now-6
   COMPONENTS
...

# Comparison
EOF
   for stride in $strides; do
     for cutoff in $cutoffs; do
       cat >>$dirname/$scaling/plumed.dat <<EOF
cv_${stride}_${cutoff//.}: LIGANDEXCHANGE ...
   GROUPA=mg
   GROUPB=op
   GROUPC=ow
   SWITCH={CUSTOM FUNC=0.5*(cos(pi*x/0.24)+1) R_0=1 D_0=0.18 D_MAX=0.42}
   VAR=_nop,_now
   FUNC=-$E*step(Y)*step(0.5-Y)*step(X-0.00001)*(1-cos(Y*pi*4))*(1-cos(X*pi))/4;X=x+$sX*(1-cos(x*pi));Y=y-$sY*(1-cos(x*pi));x=$coef*_nop-$coef*_now+6;y=$coef*_nop+$coef*_now-6
   NLIST
   NL_CUTOFF_AB=${cutoff}
   NL_CUTOFF_AC=${cutoff}
   NL_STRIDE=${stride}
   COMPONENTS
...

EOF
      done
   done
   cat >>$dirname/$scaling/plumed.dat <<EOF
scaledcv: COMBINE ARG=cv.tot COEFFICIENTS=$scaling PERIODIC=NO

uwall: UPPER_WALLS ARG=dop AT=0.8 KAPPA=50
bias: BIASVALUE ARG=scaledcv

#####################################################################
### Output ##########################################################
#####################################################################

PRINT ARG=dop,(cv\.n.*),(cv\.func.*),cv.tot,scaledcv,uwall.bias,bias.bias STRIDE=500 FILE=COLVAR
PRINT ARG=cv.tot,(cv_.*\.tot) STRIDE=1 FILE=NL_COLVAR
EOF
  done
done
