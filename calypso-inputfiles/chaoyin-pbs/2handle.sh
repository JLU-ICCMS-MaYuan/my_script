#!/bin/bash
#PBS -N  opt
#PBS -q  liuhy
#PBS -l nodes=1:ppn=28
#PBS -j oe
#PBS -V

source /public/home/mayuan/intel/oneapi/setvars.sh

ulimit -s unlimited
export I_MPI_ADJUST_REDUCE=3
export MPIR_CVAR_COLL_ALIAS_CHECK=0
export I_MPI_FABRICS=shm
export MKL_DEBUG_CPU_TYPE=5

for i in {1..300}; do
cd $i
for j in {1..4};do
cp INCAR_$j INCAR
timeout 14400s mpirun -np  28 /public/home/mayuan/software/vasp.6.1.0/bin/vasp_std > vasp.log_${j} 2>&1
wait
cp CONTCAR POSCAR
done
cd ..
done

