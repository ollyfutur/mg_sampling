#!/bin/bash

ions="MG"
replicas="rep??"

for ion in $ions; do
for rep in $replicas; do
   ion_motif_analysis $rep/input.gro $rep/run.xtc -o $rep/${ion}_OW_motifs.csv -n 8 --sel "name ${ion}" --ligsel 'nucleic or name NA MG CL OW' --inner_cutoff 2.8 --outer_cutoff 5.0 --step 10
done
done
