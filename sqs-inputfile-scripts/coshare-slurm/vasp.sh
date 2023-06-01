#!/bin/sh                           
#SBATCH  --job-name=mayqe                      
#SBATCH  --output=log.out                       
#SBATCH  --error=log.err                       
#SBATCH  --partition=normal
#SBATCH  --nodes=1                          
#SBATCH  --ntasks=64
#SBATCH  --ntasks-per-node=64
#SBATCH  --cpus-per-task=1                         

#source /work/home/may/intel/oneapi/setvars.sh --force      
#source /work/home/mayuan/intel/oneapi/setvars.sh --force      

export I_MPI_FABRICS=shm
export MKL_DEBUG_CPU_TYPE=5
source /public/env/mpi_intelmpi-2021.3.0.sh
source /public/env/compiler_intel-compiler-2021.3.0.sh
ulimit -s unlimited

mpirun -np 64 /public/software/apps/vasp/intelmpi/5.4.4/bin/vasp_std > vasp.log 2>&1