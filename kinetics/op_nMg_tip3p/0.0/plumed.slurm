#!/bin/bash
#SBATCH --job-name=kinetics
#SBATCH -N 1
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=8
#SBATCH --ntasks-per-node=1
#SBATCH --account=Sis25_bussi
#SBATCH --partition=boost_usr_prod
#SBATCH --qos=boost_qos_dbg

module purge
module load gromacs/2023/2023.5-plumed-2.10b

$GMX make_ndx -f ../equil.gro -o ../index.ndx <<EOF
ri 2 & a O1P O2P
"SOL" & a OW
q
EOF
sed -i "s/r_2_&_O1P_O2P/O_Phosphate/g" ../index.ndx &&
sed -i "s/SOL_&_OW/O_Water/g" ../index.ndx &&

 
plumed driver --ixtc run.xtc --plumed plumed_analysis.dat --trajectory-stride 0 --timestep 0.002
