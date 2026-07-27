#!/bin/bash

ion_motif_analysis 9CBU_fixed.pdb --sel "name MG" --ligsel "nucleic or name MG" --inner_cutoff 2.8 --outer_cutoff 5 --output-distances -o MG_motifs.parquet --standardize
# ion_motif_analysis rep00/equil.gro --sel "name MG" --ligsel "nucleic or name MG" --inner_cutoff 2.5 --outer_cutoff 5 --output_distances -o rep00/ion_motifs.csv
