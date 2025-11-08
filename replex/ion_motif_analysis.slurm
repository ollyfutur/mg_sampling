#!/bin/bash
#SBATCH --job-name=motifanal
#SBATCH -N 2
#SBATCH --time=24:00:00
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G
#SBATCH --ntasks-per-node=8
#SBATCH --gres=gpu:0
#SBATCH --account=Sis25_bussi
#SBATCH --partition=boost_usr_prod
#SBATCH --qos=boost_qos_dbg

ions="MG NA CL"
replicas="rep??"

for ion in $ions; do
for rep in $replicas; do
  srun -N 1 -n 1 -c 4 ion_motif_analysis $rep/input.gro $rep/run.xtc -o $rep/${ion}_motifs.csv -n 4 --sel "name ${ion}" --ligsel 'nucleic or name NA MG CL' --inner_cutoff 2.8 --outer_cutoff 5.0 &
done
done
wait
