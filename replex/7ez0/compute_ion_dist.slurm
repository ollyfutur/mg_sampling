#!/bin/bash
#SBATCH --job-name=iondist
#SBATCH -N 2
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=32
#SBATCH --gres=gpu:4
#SBATCH --account=EUHPC_R02_136
#SBATCH --partition=boost_usr_prod

module purge
module load gromacs/2023/2023.5-plumed-2.10b

ions="MG NA CL"
otypes="O_Phosphate O_Base"

replicas="rep??"
for rep in $replicas; do
  for otype in $otypes; do
    for ion in $ions; do
      srun --exclusive -N1 -n1 gmx_mpi pairdist -f $rep/run.xtc -s $rep/run.tpr -n $rep/../index.ndx -cutoff 2 -refgrouping none -selgrouping all -ref $otype -sel $ion -o $rep/dist_${otype}_${ion}.xvg &
      srun --exclusive -N1 -n1 gmx_mpi pairdist -f $rep/run.xtc -s $rep/run.tpr -n $rep/../index.ndx -cutoff 2 -refgrouping none -selgrouping all -ref $ion -sel $otype -o $rep/dist_${ion}_${otype}.xvg &
    done
  done
done
wait
