#!/bin/bash
# Job name:
#SBATCH --job-name=run_12925_correlation_vary_gamma
#
# Partition:
#SBATCH --partition=cpu-preempt
#
# Processors:
#SBATCH --ntasks=96
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
#SBATCH -o run_12925_correlation_vary_gamma.o%j
#------------------------------------------------------

echo "Starting job run_12925_correlation_vary_gamma"

module load launcher

export PYTHONNOUSERSITE=1 # ADDED TO AVOID MIXING CONDA ON CLUSTER W/ MY LOCAL PACKAGES; ADD THIS FROM NOW ON
export LAUNCHER_PLUGIN_DIR=$LAUNCHER_DIR/plugins
export LAUNCHER_RMI=SLURM
export LAUNCHER_JOB_FILE=run_12925_correlation_vary_gamma.task
export LAUNCHER_SCHED=interleaved

$LAUNCHER_DIR/paramrun

echo "Job rrun_12925_correlation_vary_gamma complete"
