#!/bin/bash
#SBATCH --job-name=iondist
#SBATCH -N 2
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=256G
#SBATCH --ntasks-per-node=32
#SBATCH --account=EUHPC_R02_136
#SBATCH --partition=boost_usr_prod
#SBATCH --qos=boost_qos_dbg

module purge
module load gromacs/2023/2023.5-plumed-2.10b

systems="7ez0 7ez0_2 7ez0_noion 7ez0_noion_2"
forcefields="mMg_tip3p mMg_opc nMg_tip3p ecc_spce"
replicas="rep00 rep05 rep11"
otypes="O_Phosphate O_Base O_Ribose N_Base"
ions="MG NA CL"

for system in $systems; do
  for forcefield in $forcefields; do
    for rep in $replicas; do
      for otype in $otypes; do
        for ion in $ions; do
          path="$system/$forcefield/$rep"
          echo "srun --exclusive -N1 -n1 -c1 gmx_mpi pairdist -f $path/run.xtc -s $path/run.tpr -n $path/../index.ndx -cutoff 2 -refgrouping none -selgrouping all -ref $otype -sel $ion -o $path/dist_${otype}_${ion}.xvg &"
          srun --exclusive -N1 -n1 -c1 --mem-per-cpu=8G gmx_mpi pairdist -f $path/run.xtc -s $path/run.tpr -n $path/../index.ndx -cutoff 2 -refgrouping none -selgrouping all -ref $otype -sel $ion -o $path/dist_${otype}_${ion}.xvg &
          echo "srun --exclusive -N1 -n1 -c1 gmx_mpi pairdist -f $path/run.xtc -s $path/run.tpr -n $path/../index.ndx -cutoff 2 -refgrouping none -selgrouping all -ref $ion -sel $otype -o $path/dist_${ion}_${otype}.xvg &"
          srun --exclusive -N1 -n1 -c1 --mem-per-cpu=8G gmx_mpi pairdist -f $path/run.xtc -s $path/run.tpr -n $path/../index.ndx -cutoff 2 -refgrouping none -selgrouping all -ref $ion -sel $otype -o $path/dist_${ion}_${otype}.xvg &
        done
      done
    done
  done
done
wait
