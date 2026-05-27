#!/bin/bash
# Job name:
#SBATCH --job-name=run_52626_krylov_norm_com
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
#SBATCH -o run_52626_krylov_norm_com.o%j
#------------------------------------------------------

echo "Starting job run_52626_krylov_norm_com"

module load launcher

export LAUNCHER_PLUGIN_DIR=$LAUNCHER_DIR/plugins
export LAUNCHER_RMI=SLURM
export LAUNCHER_JOB_FILE=run_52626_krylov_norm_com.task
export LAUNCHER_SCHED=interleaved

$LAUNCHER_DIR/paramrun

echo "Job rrun_52626_krylov_norm_com complete"
