#!/bin/bash
#SBATCH --job-name=iondist
#SBATCH --time=00:30:00
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G
#SBATCH --gres=gpu:0
#SBATCH --account=Sis25_bussi
#SBATCH --partition=boost_usr_prod
#SBATCH --qos=boost_qos_dbg

module purge
module load gromacs/2023/2023.5-plumed-2.10b

ions="MG NA"
otypes="O_Phosphate"

replicas="rep??"
for rep in $replicas; do
  for otype in $otypes; do
    for ion in $ions; do
      srun -n1 gmx_mpi pairdist -f $rep/run.xtc -s $rep/run.tpr -n $rep/../index.ndx -cutoff 2 -refgrouping none -selgrouping all -ref $otype -sel $ion -o $rep/dist_${otype}_${ion}.xvg &
      srun -n1 gmx_mpi pairdist -f $rep/run.xtc -s $rep/run.tpr -n $rep/../index.ndx -cutoff 2 -refgrouping none -selgrouping all -ref $ion -sel $otype -o $rep/dist_${ion}_${otype}.xvg &
    done
  done
done
wait
