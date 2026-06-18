#!/bin/bash
# Job name:
#SBATCH --job-name=run_6926_L4_f_sweep
#
# Partition:
#SBATCH --partition=cmt
#
# Processors:
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=128
#
# Wall clock limit:
#SBATCH --time=3-00:00:00
#
# Email at end:
#SBATCH --mail-user=sxr190099@utdallas.edu
#SBATCH --mail-type=END
#
# Specify output files:
#SBATCH -o run_6926_L4_f_sweep.o%j
#------------------------------------------------------

echo "Starting job run_6926_L4_f_sweep"

export OMP_NUM_THREADS=128
export MKL_NUM_THREADS=128
export OPENBLAS_NUM_THREADS=128
#export NUMEXPR_NUM_THREADS=1

#module load miniconda
#source ~/.bashrc
#conda activate myenv


module load launcher

export LAUNCHER_PLUGIN_DIR=$LAUNCHER_DIR/plugins
export LAUNCHER_RMI=SLURM
export LAUNCHER_JOB_FILE=run_6926_L4_f_sweep.task
export LAUNCHER_SCHED=interleaved

$LAUNCHER_DIR/paramrun

echo "Job rrun_6926_L4_f_sweep complete"
