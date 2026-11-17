#!/bin/bash
# Submit the ion counts of all systems (one array task per replicate).
cd "$PROJECT_DIR/analysis"

sbatch --array=0-1 ion_counts.slurm replex/7ez0/mMg_tip3p
sbatch --array=0-1 ion_counts.slurm replex/7ez0_noion/mMg_tip3p
sbatch --array=0-7 ion_counts.slurm plain_md/7ez0/mMg_tip3p
sbatch --array=0-7 ion_counts.slurm plain_md/7ez0/nMg_tip3p
sbatch --array=0-7 ion_counts.slurm plain_md/7ez0/nMg_tip3p_soft
sbatch --array=0-7 ion_counts.slurm plain_md/7ez0_noion/mMg_tip3p
sbatch --array=0-7 ion_counts.slurm plain_md/7ez0_noion/nMg_tip3p
sbatch --array=0-7 ion_counts.slurm plain_md/7ez0_nomg/mMg_tip3p
