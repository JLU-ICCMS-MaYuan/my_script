#!/bin/sh                   
#PBS -N    mayqe                                    
#PBS -q    liuhy         
#PBS -l    nodes=1:ppn=28               
#PBS -j    oe                                      
#PBS -V  
source /public/home/mayuan/intel/oneapi/setvars.sh --force
export I_MPI_ADJUST_REDUCE=3
export MPIR_CVAR_COLL_ALIAS_CHECK=0
ulimit -s unlimited

cd $PBS_O_WORKDIR

for i in {1..5}; do
    killall -9 vasp_std
    cp INCAR_$i INCAR
    mpirun -n 28 /public/home/mayuan/software/vasp.6.1.0/bin/vasp_std > vasp.log_$i 2>&1
    cp OUTCAR OUTCAR_$i
    cp CONTCAR POSCAR
    cp CONTCAR CONTCAR_$i
done
