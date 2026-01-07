#!/bin/bash
#SBATCH --job-name=fwddens
#SBATCH -N 1
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --ntasks-per-node=1
#SBATCH --gres=gpu:1
#SBATCH --account=Sis25_bussi
#SBATCH --partition=boost_usr_prod
#SBATCH --time=00:30:00
#SBATCH --qos=boost_qos_dbg

srun forward_density $topfile --traj $trajfile --start $start --end $end --step $step --weights $weightfile \
   --map emd_42499.mrc --sel 'not type H' --outfile $outfile --tempfactor 0 \
   --ref $reffile --align-sel 'not type H and nucleic and not (resid 22 and name P OP1 OP2)' --vox-around 6 "nucleic" \
   --processes 6 -v 
