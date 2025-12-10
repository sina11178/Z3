#!/bin/bash
# Job name:
#SBATCH --job-name=run_12925_correlation_vary_gamma
#
# Partition:
#SBATCH --partition=cmt
#
# Processors:
#SBATCH --ntasks=128
#SBATCH --nodes=1
#
# Wall clock limit:
#SBATCH --time=3-00:00:00
#
# Email at end:
#SBATCH --mail-user=sxr190099@utdallas.edu
#SBATCH --mail-type=END
#
# Specify output files:
#SBATCH -o run_12925_correlation_vary_gamma.o%j
#------------------------------------------------------

echo "Starting job run_12925_correlation_vary_gamma"

module load launcher

export LAUNCHER_PLUGIN_DIR=$LAUNCHER_DIR/plugins
export LAUNCHER_RMI=SLURM
export LAUNCHER_JOB_FILE=run_112125_autocorrelation_varying_gamma_L4_to_L10.task
export LAUNCHER_SCHED=interleaved

$LAUNCHER_DIR/paramrun

echo "Job rrun_12925_correlation_vary_gamma complete"
