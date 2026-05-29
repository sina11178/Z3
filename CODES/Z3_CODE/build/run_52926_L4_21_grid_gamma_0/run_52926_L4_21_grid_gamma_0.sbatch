#!/bin/bash
# Job name:
#SBATCH --job-name=run_52626_autocorr_vary_gamma
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
#SBATCH -o run_52626_autocorr_vary_gamma.o%j
#------------------------------------------------------

echo "Starting job run_52926_L4_21_grid_gamma_0"

#module load miniconda
#source ~/.bashrc
#conda activate myenv


module load launcher

export LAUNCHER_PLUGIN_DIR=$LAUNCHER_DIR/plugins
export LAUNCHER_RMI=SLURM
export LAUNCHER_JOB_FILE=run_52926_L4_21_grid_gamma_0.task
export LAUNCHER_SCHED=interleaved

$LAUNCHER_DIR/paramrun

echo "Job rrun_52926_L4_21_grid_gamma_0 complete"
