#!/bin/bash
# Job name:
#SBATCH --job-name=run_60126_L8_grid31_gamma_0.5_missing_pixels
#
# Partition:
#SBATCH --partition=cpu-preempt
#
# Processors:
#SBATCH --nodes=3
#SBATCH --ntasks=222
#SBATCH --cpus-per-task=2
#SBATCH --distribution=block:block
#
# Wall clock limit:
#SBATCH --time=3-00:00:00
#
# Email at end:
#SBATCH --mail-user=sxr190099@utdallas.edu
#SBATCH --mail-type=END
#
# Specify output files:
#SBATCH -o run_60126_L8_grid31_gamma_0.5_missing_pixels.o%j
#------------------------------------------------------

echo "Starting job run_60126_L8_grid31_gamma_0.5_missing_pixels"


export OMP_NUM_THREADS=2
export MKL_NUM_THREADS=2
export OPENBLAS_NUM_THREADS=2
#export NUMEXPR_NUM_THREADS=1


#module load miniconda
#source ~/.bashrc
#conda activate myenv


module load launcher

export LAUNCHER_PLUGIN_DIR=$LAUNCHER_DIR/plugins
export LAUNCHER_RMI=SLURM
export LAUNCHER_JOB_FILE=run_60126_L8_grid31_gamma_0.5_missing_pixels.task
export LAUNCHER_SCHED=interleaved

$LAUNCHER_DIR/paramrun

echo "Job rrun_60126_L8_grid31_gamma_0.5_missing_pixels complete"
