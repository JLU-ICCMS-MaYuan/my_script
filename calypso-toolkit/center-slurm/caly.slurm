#!/bin/sh
#SBATCH  --job-name=calypso
#SBATCH  --output=log.out.%j
#SBATCH  --error=log.err.%j
#SBATCH  --partition=normal
#SBATCH  --nodes=1
#SBATCH  --ntasks=64
#SBATCH  --ntasks-per-node=64
#SBATCH  --cpus-per-task=1

#source /work/env/intel2018
#conda activate py2
source /public/env/mpi_intelmpi-2021.3.0.sh
source /public/env/compiler_intel-compiler-2021.3.0.sh

ulimit -s unlimited
export I_MPI_ADJUST_REDUCE=3
export MPIR_CVAR_COLL_ALIAS_CHECK=0
export I_MPI_FABRICS=shm
export MKL_DEBUG_CPU_TYPE=5

./calypso.x > caly.log 
