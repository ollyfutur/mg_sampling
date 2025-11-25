#!/bin/bash
#SBATCH --job-name=align
#SBATCH -N 1
#SBATCH --time=24:00:00
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=32
#SBATCH --ntasks-per-node=1
#SBATCH --account=Sis25_bussi
#SBATCH --partition=boost_usr_prod
#SBATCH --qos=boost_qos_dbg

# sims="plain_md_* replex_*"
# sims="replex_7ez0_noion_mMg_tip3p_*"
# sims="plain_md_7ez0_nMg_tip3p"
sims="plain_md_7ez0_nMg_tip3p_soft plain_md_7ez0_noion_nMg_tip3p"

for sim in $sims; do
   mkdir -p $sim/aligned_residue
   srun -N1 -n1 -c32 align $sim/input.gro --mode 'residue' --traj $sim/extracted.xtc \
      --ref 9cbu.pdb \
      --sel 'not type H and nucleic' \
      --refsel 'not type H and nucleic and not (resid 22 and name P O1P O2P)' \
      --ralign 12 \
      --rout 8 \
      --reorder \
      --outsel 'not type H' \
      --output $sim/aligned_residue/aligned.pdb \
      --verbose \
      --processes 32
done
