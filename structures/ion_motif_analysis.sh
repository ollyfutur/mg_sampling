#!/bin/bash

ions="MG"
replicas="rep??"

for ion in $ions; do
for rep in $replicas; do
   # Before any minimization
   ion_motif_analysis $rep/ions2.gro -o $rep/${ion}_motifs.parquet --sel "name ${ion}" --ligsel 'nucleic or name NA MG CL' --inner_cutoff 2.8 --outer_cutoff 5.0 --standardize &
   # After equilibration
   ion_motif_analysis $rep/equil.gro -o $rep/${ion}_motifs_equil.parquet --sel "name ${ion}" --ligsel 'nucleic or name NA MG CL' --inner_cutoff 2.8 --outer_cutoff 5.0 --standardize &
done
done
wait
