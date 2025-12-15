#!/bin/bash
# Job name:
#SBATCH --job-name=run_12125_krylov_L4_to_L10_normal_comm
#
# Partition:
#SBATCH --partition=normal
#
# Processors:
#SBATCH --ntasks=64
#SBATCH --nodes=2
#
# Wall clock limit:
#SBATCH --time=3:00:00
#
# Email at end:
#SBATCH --mail-user=sxr190099@utdallas.edu
#SBATCH --mail-type=END
#
# Specify output files:
#SBATCH -o run_12125_krylov_L4_to_L10_normal_comm.o%j
#------------------------------------------------------

echo "Starting job run_12125_krylov_L4_to_L10_normal_comm"

module load launcher

export LAUNCHER_PLUGIN_DIR=$LAUNCHER_DIR/plugins
export LAUNCHER_RMI=SLURM
export LAUNCHER_JOB_FILE=run_12125_krylov_L4_to_L10_normal_comm.task
export LAUNCHER_SCHED=interleaved

$LAUNCHER_DIR/paramrun

echo "Job rrun_12125_krylov_L4_to_L10_normal_comm complete"
