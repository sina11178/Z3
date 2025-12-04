#!/bin/bash
# Job name:
#SBATCH --job-name=run_112525_autocorrelation_varying_gamma_L9
#
# Partition:
#SBATCH --partition=cmt
#
# Processors:
#SBATCH --ntasks=64
#SBATCH --nodes=1
#
# Wall clock limit:
#SBATCH --time=7-00:00:00
#
# Email at end:
#SBATCH --mail-user=sxr190099@utdallas.edu
#SBATCH --mail-type=END
#
# Specify output files:
#SBATCH -o run_112525_autocorrelation_varying_gamma_L9.o%j
#------------------------------------------------------

echo "Starting job run_112525_autocorrelation_varying_gamma_L9"

module load launcher

export LAUNCHER_PLUGIN_DIR=$LAUNCHER_DIR/plugins
export LAUNCHER_RMI=SLURM
export LAUNCHER_JOB_FILE=run_112525_autocorrelation_varying_gamma_L9.task
export LAUNCHER_SCHED=interleaved

$LAUNCHER_DIR/paramrun

echo "Job rrun_112525_autocorrelation_varying_gamma_L9 complete"
