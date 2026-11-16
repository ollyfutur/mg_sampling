#!/bin/bash
# Submit the forward density maps of all systems (one array task per replicate).
cd "$PROJECT_DIR/density"

sbatch --array=0-1 density_maps.slurm replex/7ez0/mMg_tip3p
sbatch --array=0-1 density_maps.slurm replex/7ez0_noion/mMg_tip3p
sbatch --array=0-7 density_maps.slurm plain_md/7ez0/mMg_tip3p
sbatch --array=0-7 density_maps.slurm plain_md/7ez0/nMg_tip3p
sbatch --array=0-7 density_maps.slurm plain_md/7ez0/nMg_tip3p_soft
sbatch --array=0-7 density_maps.slurm plain_md/7ez0_noion/mMg_tip3p
sbatch --array=0-7 density_maps.slurm plain_md/7ez0_noion/nMg_tip3p
sbatch --array=0-7 density_maps.slurm plain_md/7ez0_nomg/mMg_tip3p
