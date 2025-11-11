#!/bin/bash
#SBATCH --job-name=fwddens
#SBATCH -N 1
#SBATCH --time=24:00:00
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=8
#SBATCH --ntasks-per-node=1
#SBATCH --gres=gpu:1
#SBATCH --account=Sis25_bussi
#SBATCH --partition=boost_usr_prod
#SBATCH --qos=boost_qos_dbg

target=$1
srun forward_density $target/aligned_*.pdb --sel 'nucleic and not type H' --map emd_42499.map -o $target/rna.mrc --tempfactor 0 --processes 6 -v
srun forward_density $target/aligned_*.pdb --sel 'name OW' --map emd_42499.map -o $target/water.mrc --tempfactor 0 --processes 6 -v
srun forward_density $target/aligned_*.pdb --sel 'name MG NA CL' --map emd_42499.map -o $target/ions.mrc --tempfactor 0 --processes 6 -v
