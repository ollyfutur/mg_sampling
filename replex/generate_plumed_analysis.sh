#!/bin/bash

paramsfile=$1

## Parameters for the cosine envelope bias
read -r E sX sY <<< $(tail -n1 $paramsfile)

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
   COMPONENTS
...

#####################################################################
### Output ##########################################################
#####################################################################

PRINT ARG=cv.* STRIDE=5000 FILE=cv.COLVAR
"""
