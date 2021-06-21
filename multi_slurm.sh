#!/bin/bash

#SBATCH --time=72:00:00   # walltime
#SBATCH --ntasks=480   # number of processor cores (i.e. tasks)
#SBATCH --mem-per-cpu=2048M   # memory per CPU core
#SBATCH --mail-user=benp175@gmail.com   # email address
#SBATCH --mail-type=END
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=FAIL

# Set the max number of threads to use for programs using OpenMP. Should be <= ppn. Does nothing if the program doesn't use OpenMP.
export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

module load python/3.7
module load gcc/8
module load openmpi/3.1 
module load python-mpi4py/3.0
export LD_LIBRARY_PATH=/apps/gcc/8.3.0/lib64

mpiexec -n 480 python ../../../src/mm_run_multi.py

exit 0
