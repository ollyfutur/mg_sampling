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

srun forward_density $target/aligned_*.pdb --sel 'not nucleic and not type DUMMY' --map emd_42499.map -o $target/without_rna.mrc --ref 9cbu.pdb --vox-around 5 'not type H and nucleic and not resid 63-88 229-245 284-295 332-335 364-402' --tempfactor 0 --chunk-size 500000

srun forward_density $target/aligned_*.pdb --sel 'not type DUMMY' --map emd_42499.map -o $target/with_rna.mrc --ref 9cbu.pdb --vox-around 5 'not type H and nucleic and not resid 63-88 229-245 284-295 332-335 364-402' --tempfactor 0 --chunk-size 500000
