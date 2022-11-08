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

source activate cage
#dftbenchmark -i ecutjobs.ini -e > test.log 2>&1 
dftbenchmark -i kspacing.ini -e > test.log 2>&1 
