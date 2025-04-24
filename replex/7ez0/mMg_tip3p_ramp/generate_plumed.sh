#!/bin/bash

scalingsfile=$1
ramp_scalingsfile=$2

# Get scaling factors
readarray -t scalings < $scalingsfile
readarray -t ramp_scalings < $ramp_scalingsfile

### Generate Plumed files
## Parameters for the cosine envelope bias
## (From optimization)
read -r E sX sY <<< $(tail -n1 bias_params.wsv)

# Needed to rescale coordination number to integer values
coef=1.039661

# Generate Plumed input file for replex
for i in "${!scalings[@]}"; do
  a="${ramp_scalings[i]}"
  lambda="${scalings[i]}"
  printf "lambda: %s / ramp lambda: %s\n" $lambda $a
  repdir=$(printf "rep%02d" $i)
  echo """
LOAD FILE=$PROJECT_DIR/ligandexchange/LigandExchange.so

rmsd: RMSD REFERENCE=../ref_subset.pdb TYPE=SIMPLE

mg: GROUP NDX_FILE=../index.ndx NDX_GROUP=MG
op: GROUP NDX_FILE=../index.ndx NDX_GROUP=O_Phosphate
ob: GROUP NDX_FILE=../index.ndx NDX_GROUP=O_Base
ox: GROUP ATOMS=op,ob
ow: GROUP NDX_FILE=../index.ndx NDX_GROUP=O_Water

cv: LIGANDEXCHANGE ...
   GROUPA=mg
   GROUPB=ox
   GROUPC=ow
   SWITCH={COSINUS R_0=0.24 D_0=0.18}
   VAR=_nox,_now
   FUNC=-$E*$lambda*step(Y)*step(0.5-Y)*step(X-0.00001)*(1-cos(Y*pi*4))*(1-cos(X*pi))/4+$a*$coef*_nox;X=x+$sX*(1-cos(x*pi));Y=y-$sY*(1-cos(x*pi));x=$coef*_nox-$coef*_now+6;y=$coef*_nox+$coef*_now-6
   NLIST
   NL_CUTOFF_AB=1.0
   NL_CUTOFF_AC=1.0
   NL_STRIDE=500
   NESTED
...

#####################################################################
### Biases ##########################################################
#####################################################################

restr: RESTRAINT ARG=rmsd AT=0 KAPPA=500

bias: BIASVALUE ARG=cv

#####################################################################
### Output ##########################################################
#####################################################################

PRINT ARG=rmsd,restr.bias STRIDE=5000 FILE=rmsd.COLVAR
PRINT ARG=cv,bias.bias STRIDE=5000 FILE=bias.COLVAR
""" > $repdir/plumed_ramp.02.dat
done
