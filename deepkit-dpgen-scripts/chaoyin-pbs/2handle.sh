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

cd $PBS_O_WORKDIR
killall -9 vasp_std


for i in {1..500}; do
    mkdir $i
    cd $i
    cp ../INCAR .
    cp ../POTCAR .
    cp ../POSCAR_$i POSCAR_$i
    cp POSCAR_$i POSCAR
    mpirun -np  28 /public/home/mayuan/software/vasp.6.1.0/bin/vasp_std > vasp.log  2>&1
    cd ..
done

