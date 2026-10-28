#!/bin/bash

paramsfile=$1
scalingsfile=$2
alphasfile=$3
# Index to take in scalingsfile and alphasfile
i=$4

## Parameters for the cosine envelope bias
read -r E sX sY <<< $(tail -n1 $paramsfile)

# Get scaling factors
readarray -t scalings < $scalingsfile
scaling=$(printf "%.6f" "${scalings[i]}")

# Get ramp slopes
readarray -t alphas < $alphasfile
alpha=$(printf "%.6f" "${alphas[i]}")

### Generate Plumed file
# Needed to rescale coordination number to integer values
coef=1.039661

echo """
LOAD FILE=$PROJECT_DIR/ligandexchange/LigandExchange.so

mg: GROUP NDX_FILE=../index.ndx NDX_GROUP=MG
op: GROUP NDX_FILE=../index.ndx NDX_GROUP=O_Phosphate
ow: GROUP NDX_FILE=../index.ndx NDX_GROUP=O_Water

# (unscaled) barrier-flattening bias
cv: LIGANDEXCHANGE ...
   GROUPA=mg
   GROUPB=op
   GROUPC=ow
   SWITCH={CUSTOM FUNC=0.5*(cos(pi*x/0.24)+1) R_0=1 D_0=0.18 D_MAX=0.42}
   VAR=_nop,_now
   FUNC=-$E*step(Y)*step(0.5-Y)*step(X-0.00001)*(1-cos(Y*pi*4))*(1-cos(X*pi))/4;X=x+$sX*(1-cos(x*pi));Y=y-$sY*(1-cos(x*pi));x=$coef*_nop-$coef*_now+6;y=$coef*_nop+$coef*_now-6
   NLIST
   NL_CUTOFF_AB=1.0
   NL_CUTOFF_AC=1.0
   NL_STRIDE=250
   NESTED
...

# redundant but needed for ramp correction (LIGANDEXCHANGE does not expose it)
_nop: COORDINATION ...
   GROUPA=mg
   GROUPB=op
   SWITCH={CUSTOM FUNC=0.5*(cos(pi*x/0.24)+1) R_0=1 D_0=0.18 D_MAX=0.42}
   NLIST
   NL_CUTOFF=1.0
   NL_STRIDE=250
...
nop: COMBINE ARG=_nop COEFFICIENTS=1.039661 PERIODIC=NO
cvtot: COMBINE ARG=cv,nop COEFFICIENTS=$scaling,$alpha PERIODIC=NO

#####################################################################
### Biases ##########################################################
#####################################################################

bias: BIASVALUE ARG=cvtot

#####################################################################
### Output ##########################################################
#####################################################################

PRINT ARG=cv,nop,bias.bias STRIDE=5000 FILE=bias.COLVAR
"""
