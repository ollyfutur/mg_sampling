#!/bin/bash
#SBATCH --job-name=trjconv
#SBATCH -N 1
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --account=Sis25_bussi
#SBATCH --partition=boost_usr_prod
#SBATCH --qos=boost_qos_dbg

$GMX check -f run.xtc > check_xtc.log
#$GMX trjconv -f run.xtc -s run.tpr -o run_light.xtc -dt 10 << EOM
#System
#EOM
