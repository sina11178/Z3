#!/bin/bash
# Job name:
#SBATCH --job-name=run_53126_L8_31_grid_gamma_1
#
# Partition:
#SBATCH --partition=cmt
#
# Processors:
#SBATCH --ntasks=256
#SBATCH --nodes=2
#
# Wall clock limit:
#SBATCH --time=3-00:00:00
#
# Email at end:
#SBATCH --mail-user=sxr190099@utdallas.edu
#SBATCH --mail-type=END
#
# Specify output files:
#SBATCH -o run_53126_L8_31_grid_gamma_1.o%j
#------------------------------------------------------

echo "Starting job run_53126_L8_31_grid_gamma_1"


export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
#export NUMEXPR_NUM_THREADS=1


#module load miniconda
#source ~/.bashrc
#conda activate myenv


module load launcher

export LAUNCHER_PLUGIN_DIR=$LAUNCHER_DIR/plugins
export LAUNCHER_RMI=SLURM
export LAUNCHER_JOB_FILE=run_53126_L8_31_grid_gamma_1.task
export LAUNCHER_SCHED=interleaved

$LAUNCHER_DIR/paramrun

echo "Job rrun_53126_L8_31_grid_gamma_1 complete"
