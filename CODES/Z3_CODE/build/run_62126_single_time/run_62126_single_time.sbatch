#!/bin/bash
# Job name:
#SBATCH --job-name=run_62126_single_time
#
# Partition:
#SBATCH --partition=cmt
#
# Processors:
#SBATCH --ntasks=3
#SBATCH --nodes=2
#SBATCH --cpus-per-task=64
#
# Wall clock limit:
#SBATCH --time=3-00:00:00
#
# Email at end:
#SBATCH --mail-user=sxr190099@utdallas.edu
#SBATCH --mail-type=END
#
# Specify output files:
#SBATCH -o run_62126_single_time.o%j
#------------------------------------------------------

echo "Starting job run_62126_single_time"

export OMP_NUM_THREADS=85
export MKL_NUM_THREADS=85
export OPENBLAS_NUM_THREADS=85
#export NUMEXPR_NUM_THREADS=1

#module load miniconda
#source ~/.bashrc
#conda activate myenv


module load launcher

export LAUNCHER_PLUGIN_DIR=$LAUNCHER_DIR/plugins
export LAUNCHER_RMI=SLURM
export LAUNCHER_JOB_FILE=run_62126_single_time.task
export LAUNCHER_SCHED=interleaved

$LAUNCHER_DIR/paramrun

echo "Job rrun_62126_single_time complete"
