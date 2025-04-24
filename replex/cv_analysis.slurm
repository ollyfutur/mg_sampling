#!/bin/bash
#SBATCH --job-name=cvanal
#SBATCH -N 1
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=24
#SBATCH --account=EUHPC_R02_136
#SBATCH --partition=boost_usr_prod

module purge
module load gromacs/2023/2023.5-plumed-2.10b

systems="7ez0 7ez0_noion"
forcefields="mMg_opc mMg_tip3p nMg_tip3p ecc_spce"
repids="00 06 11"

workdir=$(pwd)

for system in $systems; do
  for forcefield in $forcefields; do
    for repid in $repids; do
      simdir="$system/$forcefield/rep$repid"
      cd $simdir
      srun -n 1 -c 1 plumed driver --plumed ../plumed.analysis.dat --ixtc run.xtc --trajectory-stride 0 &
      cd $workdir
    done
  done
done
wait
