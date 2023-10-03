#!/bin/sh                                     
#SBATCH  --job-name=calypso
#SBATCH  --output=log.calypso.out.%j             
#SBATCH  --error=log.calypso.err.%j              
#SBATCH  --partition=normal                       
#SBATCH  --nodes=1                            
#SBATCH  --ntasks=64
#SBATCH  --ntasks-per-node=64             
#SBATCH  --cpus-per-task=1                    

source /public/env/mpi_intelmpi-2021.3.0.sh
source /public/env/compiler_intel-compiler-2021.3.0.sh

ulimit -s unlimited
export I_MPI_ADJUST_REDUCE=3
export MPIR_CVAR_COLL_ALIAS_CHECK=0
export I_MPI_FABRICS=shm
export MKL_DEBUG_CPU_TYPE=5

for i in {1..500}; do
    mkdir $i
    cd $i
    cp ../INCAR .
    cp ../POTCAR .
    cp ../POSCAR_$i POSCAR_$i
    cp POSCAR_$i POSCAR
    mpirun -n 64 /public/software/apps/vasp/intelmpi/5.4.4/bin/vasp_std > vasp.log 2>&1
    cd ..
done

