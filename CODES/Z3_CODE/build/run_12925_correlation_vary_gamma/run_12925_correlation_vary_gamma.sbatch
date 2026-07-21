#!/bin/bash
# Job name:
#SBATCH --job-name=run_12925_correlation_vary_gamma
#
# Partition:
#SBATCH --partition=cmt
#
# Processors:
#SBATCH --ntasks=60
#SBATCH --nodes=2
#SBATCH --cpus-per-task=4
#
# Wall clock limit:
#SBATCH --time=2-00:00:00
#
# Email at end:
#SBATCH --mail-user=sxr190099@utdallas.edu
#SBATCH --mail-type=END
#
# Specify output files:
#SBATCH -o run_12925_correlation_vary_gamma.o%j
#------------------------------------------------------

echo "Starting job run_12925_correlation_vary_gamma"


export OMP_NUM_THREADS=4
export MKL_NUM_THREADS=4
export OPENBLAS_NUM_THREADS=4
#export NUMEXPR_NUM_THREADS=1

module load launcher

export LAUNCHER_PLUGIN_DIR=$LAUNCHER_DIR/plugins
export LAUNCHER_RMI=SLURM
export LAUNCHER_JOB_FILE=run_12925_correlation_vary_gamma.task
export LAUNCHER_SCHED=interleaved

$LAUNCHER_DIR/paramrun

echo "Job rrun_12925_correlation_vary_gamma complete"
