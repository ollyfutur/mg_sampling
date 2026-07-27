#!/bin/bash

ions="MG NA"
replicas="rep??"

for ion in $ions; do
for rep in $replicas; do
   if [ ! -f $rep/${ion}_motifs.parquet ]; then
      ion_motif_analysis $rep/input.gro $rep/run.xtc -o $rep/${ion}_motifs.parquet --sel "name ${ion}" --ligsel 'nucleic or name NA MG CL' --inner_cutoff 2.8 --outer_cutoff 5.0 --standardize &
   fi
done
done
wait
