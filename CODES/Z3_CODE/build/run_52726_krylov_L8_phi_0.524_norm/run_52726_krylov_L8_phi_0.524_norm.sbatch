#!/bin/bash
# Job name:
#SBATCH --job-name=run_52726_krylov_L8_phi_0.524_norm
#
# Partition:
#SBATCH --partition=normal
#
# Processors:
#SBATCH --ntasks=32
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
#SBATCH -o run_52726_krylov_L8_phi_0.524_norm.o%j
#------------------------------------------------------

echo "Starting job run_52726_krylov_L8_phi_0.524_norm"

module load launcher

export LAUNCHER_PLUGIN_DIR=$LAUNCHER_DIR/plugins
export LAUNCHER_RMI=SLURM
export LAUNCHER_JOB_FILE=run_52726_krylov_L8_phi_0.524_norm.task
export LAUNCHER_SCHED=interleaved

$LAUNCHER_DIR/paramrun

echo "Job rrun_52726_krylov_L8_phi_0.524_norm complete"
