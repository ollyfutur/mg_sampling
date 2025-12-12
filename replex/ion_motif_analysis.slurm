#!/bin/bash
#SBATCH --job-name=motif
#SBATCH -N 1
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --ntasks-per-node=8
#SBATCH --account=Sis25_bussi
#SBATCH --partition=boost_usr_prod

ions="MG NA"
sims="7ez0/mMg_tip3p/replicate0 7ez0/mMg_tip3p/replicate1 7ez0_noion/mMg_tip3p/replicate0 7ez0_noion/mMg_tip3p/replicate1"

for ion in $ions; do
   for sim in $sims; do
      replicas="$sim/rep??"
      for rep in $replicas; do
         srun -n1 -c4 ion_motif_analysis $rep/input.gro $rep/run.xtc -o $rep/${ion}_motifs.parquet --sel "name ${ion}" --ligsel 'nucleic or name NA MG CL' --inner_cutoff 2.8 --outer_cutoff 0 -n 4 &
      done
      wait
   done
done
