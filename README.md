This is the repository containing simulation inputs and analysis code for (Enhanced sampling and cryo-EM data resolve magnesium binding to RNA)[https://doi.org/10.48550/arXiv.2607.07186].

# Structure of the repository


# Use the accelerating bias

1) Requirements:
- GROMACS patched with PLUMED >=2.10
- The following modules:
```bash
git clone https://github.com/ollyfutur/ligandexchange.git
git clone https://github.com/ollyfutur/replex_biasvalue.git
```
- Compile the shared objects (make sure you have PLUMED properly loaded in your environment)
```bash
cd ligandexchange
make
cd ..
cd replex_biasvalue
make
cd ..
```

2) Prepare your system:
- The bias was parameterized to be used with parmBSC0χOL3, TIP3P and microMg
- You may use `forcefields/amber14sb_mMg_tip3p.ff`
> [!WARNING]
> If you keep Mg atoms from the input (e.g. PDB) structure, make sure the topology file uses the `mMg` atom type:
> ```bash
> sed -i 's/\(^\s*[0-9]*\s*\)MG /\1mMg/g' topol_Ion_chain_N2.itp
> ```

3) Prepare the bias parameters:


4) Prepare PLUMED input file using the following template

```plumed
#TODO: change path to 
LOAD FILE={path/to/ligandexchange/LigandExchange.so}
LOAD FILE={path/to/replex_biasvalue/ReplexBiasValue.so}

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
   FUNC=-38.307903*step(Y)*step(0.5-Y)*step(X-0.00001)*(1-cos(Y*pi*4))*(1-cos(X*pi))/4;X=x+0.10098*(1-cos(x*pi));Y=y-0.287523*(1-cos(x*pi));x=1.039661*_nop-1.039661*_now+6;y=1.039661*_nop+1.039661*_now-6
   NLIST
   NL_CUTOFF_AB=1.0
   NL_CUTOFF_AC=1.0
   NL_STRIDE=250
...

# ramp correction
_nop: COORDINATION ...
   GROUPA=mg
   GROUPB=op
   SWITCH={CUSTOM FUNC=0.5*(cos(pi*x/0.24)+1) R_0=1 D_0=0.18 D_MAX=0.42}
   NLIST
   NL_CUTOFF=1.0
   NL_STRIDE=250
...
nop: COMBINE ARG=_nop COEFFICIENTS=1.039661 PERIODIC=NO

#####################################################################
### Biases ##########################################################
#####################################################################

bias: REPLEX_BIASVALUE ...
   ARG=cv,nop
   COEFFICIENTS={
#TODO: insert coefficients here
   }
   PACE=500
   NEX=100
...

#####################################################################
### Output ##########################################################
#####################################################################

PRINT ARG=cv,nop,bias.bias,bias.index STRIDE=5000 FILE=bias.COLVAR

```

# How to reproduce the results

