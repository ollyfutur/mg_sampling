#!/bin/bash
#SBATCH --job-name=motifanal
#SBATCH -N 1
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=16
#SBATCH --ntasks-per-node=1
#SBATCH --account=EUHPC_R02_136
#SBATCH --partition=boost_usr_prod

module purge

simdirs="*/*/rep??"

for simdir in $simdirs; do
  srun ion_motif_analysis $simdir/input.gro $simdir/run.xtc -o $simdir/ion_analysis.csv -n $SLURM_CPUS_PER_TASK --step 10
done
