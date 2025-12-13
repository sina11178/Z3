#!/bin/bash
# Job name:
#SBATCH --job-name=run_121325_L8_gam0.4_phi_0.524
#
# Partition:
#SBATCH --partition=gpu-preempt
#
# Processors:
#SBATCH --ntasks=16
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
#SBATCH -o run_121325_L8_gam0.4_phi_0.524.o%j
#------------------------------------------------------

echo "Starting job run_121325_L8_gam0.4_phi_0.524"

module load launcher

export LAUNCHER_PLUGIN_DIR=$LAUNCHER_DIR/plugins
export LAUNCHER_RMI=SLURM
export LAUNCHER_JOB_FILE=run_121325_L8_gam0.4_phi_0.524.task
export LAUNCHER_SCHED=interleaved

$LAUNCHER_DIR/paramrun

echo "Job rrun_121325_L8_gam0.4_phi_0.524 complete"
