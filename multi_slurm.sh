#!/bin/bash

#SBATCH --time=20:00:00   # walltime
#SBATCH --ntasks=40   # number of processor cores (i.e. tasks)
# #SBATCH --nodes=10   # number of nodes
#SBATCH --mem-per-cpu=100M   # memory per CPU core
#SBATCH --mail-user=dallinspencer@gmail.com   # email address
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL


# Set the max number of threads to use for programs using OpenMP. Should be <= ppn. Does nothing if the program doesn't use OpenMP.
export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

#module load miniconda3/4.6
#source activate mmenv
module load python/3.7
module load gcc/8

#module load llvm/7
module load openmpi/3.1 
module load python-mpi4py/3.0
export LD_LIBRARY_PATH=/apps/gcc/8.3.0/lib64

mpiexec -n 40 python ../../../src/mm_run_multi.py

# python3 ../../../src/mm_run.py

exit 0
