#!/bin/bash
#SBATCH --job-name=iondist
#SBATCH -N 1
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=24
#SBATCH --gres=gpu:3
#SBATCH --account=EUHPC_R02_136
#SBATCH --partition=boost_usr_prod
#SBATCH --qos=boost_qos_dbg

module purge
module load gromacs/2023/2023.5-plumed-2.10b

ions="MG NA CL"
otypes="O_Phosphate O_Base O_Ribose N_Base"

replicas="rep00 rep05 rep11"
for rep in $replicas; do
  for otype in $otypes; do
    for ion in $ions; do
      srun --exclusive -N1 -n1 gmx_mpi pairdist -f $rep/run.xtc -s $rep/run.tpr -n $rep/../index.ndx -cutoff 2 -refgrouping none -selgrouping all -ref $otype -sel $ion -o $rep/dist_${otype}_${ion}.xvg &
      srun --exclusive -N1 -n1 gmx_mpi pairdist -f $rep/run.xtc -s $rep/run.tpr -n $rep/../index.ndx -cutoff 2 -refgrouping none -selgrouping all -ref $ion -sel $otype -o $rep/dist_${ion}_${otype}.xvg &
    done
  done
  wait
done
