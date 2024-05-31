#!/bin/sh                                     
#SBATCH  --job-name=vasp
#SBATCH  --output=log.vasp.out.%j             
#SBATCH  --error=log.vasp.err.%j              
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


for i in {1..4}; do 
cp INCAR_$i INCAR
timeout 14400s mpirun -n 64 /public/software/apps/vasp/intelmpi/5.4.4/bin/vasp_std > vasp.log_${i} 2>&1
wait
cp CONTCAR POSCAR
done

rm -rf CHG* WAVECAR

