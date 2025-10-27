#!/bin/bash
#SBATCH --job-name=align
#SBATCH -N 1
#SBATCH --time=24:00:00
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=256G
#SBATCH --ntasks-per-node=1
#SBATCH --account=Sis25_bussi
#SBATCH --partition=boost_usr_prod
#SBATCH --qos=boost_qos_dbg

mkdir -p aligned_global
srun -N1 -n1 -c32 --mem-per-cpu=8G align input.gro --mode 'global' --traj extracted.xtc \
   --ref ../9cbu.pdb \
   --sel 'not type H and nucleic and not resid 63-88 229-245 284-295 332-335 364-402' \
   --refsel 'not type H and nucleic and not resid 63-88 229-245 284-295 332-335 364-402 and not (resid 22 and name P O1P O2P)' \
   --ralign 12 \
   --rout 8 \
   --reorder \
   --outsel 'not type H' \
   --output aligned_global/aligned.pdb \
   -p 32
