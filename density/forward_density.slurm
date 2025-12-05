#!/bin/bash
#SBATCH --job-name=fwddens
#SBATCH -N 1
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --ntasks-per-node=1
#SBATCH --gres=gpu:1
#SBATCH --account=Sis25_bussi
#SBATCH --partition=boost_usr_prod

for target in plain*/aligned_residue replex*/aligned_residue; do
   files="$target/aligned_?.pdb $target/aligned_??.pdb $target/aligned_???.pdb $target/aligned_????.pdb"
   #files="$target/aligned_?.pdb" 
   srun forward_density $files --sel 'all' --map emd_42499.mrc -o $target/all.mrc --tempfactor 0 --processes 6 -v --frame-out --vox-around 5 "nucleic" --ref 9CBU.pdb --vox-indices-out $target/vox_indices.npy
done
# srun forward_density $target/aligned_*.pdb --sel 'nucleic' --map emd_42499.mrc -o $target/rna.mrc --tempfactor 0 --processes 6 -v --chunk-size 32768
# srun forward_density $target/aligned_*.pdb --sel 'name MG' --map emd_42499.mrc -o $target/mg.mrc --tempfactor 0 --processes 6 -v
# srun forward_density $target/aligned_*.pdb --sel 'name NA' --map emd_42499.mrc -o $target/na.mrc --tempfactor 0 --processes 6 -v
