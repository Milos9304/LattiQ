#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1      ##This should be scheduled automatically
#SBATCH --cpus-per-task=48       ##This should be scheduled automatically
#SBATCH --exclusive

#SBATCH --nvram-options=none

# set number of qubits
#SBATCH --job-name=VQESuperInstance
NUM_QUBITS=30

# OMP settings
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PROC_BIND=TRUE
export OMP_PLACES=cores

# Report run information
echo "SLURM_JOB_NUM_NODES=$SLURM_JOB_NUM_NODES"
echo "SLURM_TASKS_PER_NODE=$SLURM_TASKS_PER_NODE"
echo "SLURM_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK"
echo "OMP_NUM_THREADS=$OMP_NUM_THREADS"

srun --distribution=block:block ../bin/fastVQA --vqe -x 0.175 -r 30 -d qary_25_50 --seed 1997
# Report job time
echo ""
echo "ROUGH WALLTIME: (mm:ss)"
squeue -h -j $SLURM_JOBID -o "%M"