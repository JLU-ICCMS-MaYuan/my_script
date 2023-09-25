#!/bin/sh                                     
#SBATCH  --job-name=calypso
#SBATCH  --output=log.calypso.out.%j             
#SBATCH  --error=log.calypso.err.%j              
#SBATCH  --partition=lhy
#SBATCH  --nodes=1                            
#SBATCH  --ntasks=48
#SBATCH  --ntasks-per-node=48
#SBATCH  --cpus-per-task=1                    

source /work/home/may/intel/oneapi/setvars.sh --force

ulimit -s unlimited
export I_MPI_ADJUST_REDUCE=3
export MPIR_CVAR_COLL_ALIAS_CHECK=0
export I_MPI_FABRICS=shm
export MKL_DEBUG_CPU_TYPE=5

for i in {1..300}; do
cd $i
for j in {1..4};do
cp INCAR_$j INCAR
timeout 14400s mpirun -np 48 /work/software/vasp.6.1.0/vasp_std > vasp.log_${j} 2>&1
wait
cp CONTCAR POSCAR
done
cd ..
done

