#!/bin/bash
#SBATCH --job-name=motifanal
#SBATCH -N 1
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=16
#SBATCH --ntasks-per-node=2
#SBATCH --account=EUHPC_R02_136
#SBATCH --partition=boost_usr_prod

module purge

systems="7ez0 7ez0_noion"
forcefields="ecc_spce"
repids="11"

for system in $systems; do
  for forcefield in $forcefields; do
    for repid in $repids; do
      simdir="$system/$forcefield/rep$repid"
      srun -n 1 -c 16 ion_motif_analysis $simdir/input.gro $simdir/run.xtc -o $simdir/ion_analysis.csv -n 16 &
    done
  done
done
wait
