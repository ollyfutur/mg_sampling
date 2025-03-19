#!/bin/bash

simdirs="$@"
ions="MG NA CL"
otypes="O_Phosphate O_Base"
for simdir in $simdirs; do
  for otype in $otypes; do
    for ion in $ions; do
      gmx_mpi pairdist -f $simdir/run.xtc -s $simdir/run.tpr -n $simdir/../index.ndx -cutoff 2 -refgrouping none -selgrouping all -ref $otype -sel $ion -o $simdir/dist_${otype}_${ion}.xvg &
      gmx_mpi pairdist -f $simdir/run.xtc -s $simdir/run.tpr -n $simdir/../index.ndx -cutoff 2 -refgrouping none -selgrouping all -ref $ion -sel $otype -o $simdir/dist_${ion}_${otype}.xvg &
    done
  done
  wait
done
