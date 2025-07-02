#!/bin/bash
# Job name:
#SBATCH --job-name=run_063025_krylov_matrices_zero_modes_L4_to_L6
#
# Partition:
#SBATCH --partition=normal
#
# Processors:
#SBATCH --ntasks=144
#SBATCH --nodes=9
#
# Wall clock limit:
#SBATCH --time=2-00:00:00
#
# Email at end:
#SBATCH --mail-user=sxr190099@utdallas.edu
#SBATCH --mail-type=END
#
# Specify output files:
#SBATCH -o run_063025_krylov_matrices_zero_modes_L4_to_L6.o%j
#------------------------------------------------------

echo "Starting job run_063025_krylov_matrices_zero_modes_L4_to_L6"

module load launcher

export LAUNCHER_PLUGIN_DIR=$LAUNCHER_DIR/plugins
export LAUNCHER_RMI=SLURM
export LAUNCHER_JOB_FILE=run_063025_krylov_matrices_zero_modes_L4_to_L6.task
export LAUNCHER_SCHED=interleaved

$LAUNCHER_DIR/paramrun

echo "Job run_063025_krylov_matrices_zero_modes_L4_to_L6 complete"
