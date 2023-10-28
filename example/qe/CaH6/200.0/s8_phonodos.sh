#!/bin/sh   
source /public/home/mayuan/intel/oneapi/setvars.sh --force
ulimit -s unlimited
mpirun -np 1 /public/home/mayuan/software/qe-7.1/bin/matdyn.x  <phonodos.in> phonodos.out  
