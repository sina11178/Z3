#!/bin/bash
# Job name:
#SBATCH --job-name=run_210807_fPhi_PhaseDiag
#
# Partition:
#SBATCH --partition=small
#
# Account:
#SBATCH -A inverted_mbl # which account to charge the SUs to
#
# Processors:
#SBATCH --ntasks=19
#SBATCH --nodes=2
#SBATCH --cpus-per-task=5
#
# Wall clock limit:
#SBATCH --time=2-00:00:00
#
# Email at end:
#SBATCH --mail-user=sxr176930@utdallas.edu
#SBATCH --mail-type=END
#
# Specify output files:
#SBATCH -o run_210807_fPhi_PhaseDiag.o%j
#------------------------------------------------------

echo "Starting job run_210807_fPhi_PhaseDiag"
echo " l==9 and f==0.1 and and 2.0<phi<4.0"
module load launcher

export LAUNCHER_PLUGIN_DIR=$LAUNCHER_DIR/plugins
export LAUNCHER_RMI=SLURM
export LAUNCHER_JOB_FILE=run_210807_fPhi_PhaseDiag.task
export LAUNCHER_SCHED=interleaved

$LAUNCHER_DIR/paramrun

echo "Job run_210807_fPhi_PhaseDiag complete"
