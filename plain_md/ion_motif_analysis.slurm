#!/bin/bash
#SBATCH --job-name=motif
#SBATCH -N 1
#SBATCH --time=24:00:00
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=4
#SBATCH --ntasks-per-node=8
#SBATCH --account=Sis25_bussi
#SBATCH --partition=boost_usr_prod
#SBATCH --qos=boost_qos_dbg

ions="MG"
replicas="rep??"

for ion in $ions; do
   for rep in $replicas; do
      srun -n1 -c4 ion_motif_analysis $rep/input.gro $rep/run.xtc -o $rep/${ion}_OW_motifs.parquet --sel "name ${ion}" --ligsel 'nucleic or name NA MG CL OW' --inner_cutoff 2.8 --outer_cutoff 5.0 --step 10 -n 4 &
   done
done
wait
