#!/bin/bash
# Job name:
#SBATCH --job-name=run_042225_autocorr_vary_gamma
#
# Partition:
#SBATCH --partition=normal
#
# Processors:
#SBATCH --ntasks=32
#SBATCH --nodes=2
#
# Wall clock limit:
#SBATCH --time=2-00:00:00
#
# Email at end:
#SBATCH --mail-user=sxr190099@utdallas.edu
#SBATCH --mail-type=END
#
# Specify output files:
#SBATCH -o run_042225_autocorr_vary_gamma.o%j
#------------------------------------------------------

echo "Starting job run_042525_autocorr_vary_gamma"

module load launcher

export LAUNCHER_PLUGIN_DIR=$LAUNCHER_DIR/plugins
export LAUNCHER_RMI=SLURM
export LAUNCHER_JOB_FILE=run_042525_autocorr_vary_gamma.task
export LAUNCHER_SCHED=interleaved

$LAUNCHER_DIR/paramrun

echo "Job rrun_042525_autocorr_vary_gamma complete"
