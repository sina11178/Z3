#!/bin/bash
# Job name:
#SBATCH --job-name=run_121125_L8_only_twisted
#
# Partition:
#SBATCH --partition=cmt
#
# Processors:
#SBATCH --ntasks=64
#SBATCH --nodes=1
#
# Wall clock limit:
#SBATCH --time=3:00:00
#
# Email at end:
#SBATCH --mail-user=sxr190099@utdallas.edu
#SBATCH --mail-type=END
#
# Specify output files:
#SBATCH -o run_121125_L8_only_twisted.o%j
#------------------------------------------------------

echo "Starting job run_121125_L8_only_twisted"

module load launcher

export LAUNCHER_PLUGIN_DIR=$LAUNCHER_DIR/plugins
export LAUNCHER_RMI=SLURM
export LAUNCHER_JOB_FILE=run_121125_L8_only_twisted.task
export LAUNCHER_SCHED=interleaved

$LAUNCHER_DIR/paramrun

echo "Job rrun_121125_L8_only_twisted complete"
