#!/bin/bash

#SBATCH --time=72:00:00   # walltime
#SBATCH --ntasks=100   # number of processor cores (i.e. tasks)
#SBATCH --nodes=4   # number of nodes
#SBATCH --mem-per-cpu=512M   # memory per CPU core
#SBATCH --mail-user=dallinspencer@gmail.com   # email address
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL


# Set the max number of threads to use for programs using OpenMP. Should be <= ppn. Does nothing if the program doesn't use OpenMP.
export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
module load python/3.8
module load mpi

mpiexec -n 100 python3 ../../../src/mm_run.py

exit 0
