#!/bin/bash
# Job name:
#SBATCH --job-name=run_52626_autocorr_vary_gamma
#
# Partition:
#SBATCH --partition=cpu-preempt
#
# Processors:
#SBATCH --ntasks=192
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


export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
#export NUMEXPR_NUM_THREADS=1


echo "Starting job optimize_calc_time_53026"

#module load miniconda
#source ~/.bashrc
#conda activate myenv


module load launcher

export LAUNCHER_PLUGIN_DIR=$LAUNCHER_DIR/plugins
export LAUNCHER_RMI=SLURM
export LAUNCHER_JOB_FILE=optimize_calc_time_53026.task
export LAUNCHER_SCHED=interleaved

$LAUNCHER_DIR/paramrun

echo "Job roptimize_calc_time_53026 complete"
