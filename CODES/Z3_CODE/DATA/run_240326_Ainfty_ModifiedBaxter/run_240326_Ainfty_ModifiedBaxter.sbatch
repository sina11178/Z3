#!/bin/bash
# Job name:
#SBATCH --job-name=run_240326_Ainfty_ModifiedBaxter
#
# Partition:
#SBATCH --partition=normal
#
# Account:
#SBATCH -A PHY23013 # which account to charge the SUs to
#
# Processors:
#SBATCH --ntasks=600
#SBATCH --nodes=5
#SBATCH --cpus-per-task=1
#
# Wall clock limit:
#SBATCH --time=2-00:00:00
#
# Email at end:
#SBATCH --mail-user=sxr176930@utdallas.edu
#SBATCH --mail-type=END
#
# Specify output files:
#SBATCH -o run_240326_Ainfty_ModifiedBaxter.o%j
#------------------------------------------------------

echo "Starting job run_240326_Ainfty_ModifiedBaxter"
echo "n2 = 50 and n3 = 1600 "
module load launcher

export LAUNCHER_PLUGIN_DIR=$LAUNCHER_DIR/plugins
export LAUNCHER_RMI=SLURM
export LAUNCHER_JOB_FILE=run_240326_Ainfty_ModifiedBaxter.task
export LAUNCHER_SCHED=interleaved

$LAUNCHER_DIR/paramrun

echo "Job run_240326_Ainfty_ModifiedBaxter complete"
