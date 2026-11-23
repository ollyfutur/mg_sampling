#!/bin/bash
# Mg²⁺ binding motifs of the deposited structure (reference sites), as for 9cbu.

ion_motif_analysis ../7EZ0_fixed.pdb --sel "name MG" --ligsel "nucleic or name MG" --inner_cutoff 2.8 --outer_cutoff 5 --output-distances -o MG_motifs.parquet --standardize
