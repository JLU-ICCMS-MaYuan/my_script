#!/bin/sh                                                                
#PBS -N    relax                                                         
#PBS -q    liuhy                                                         
#PBS -l    nodes=1:ppn=28                                                
#PBS -j    oe                                                            
#PBS -V                                                                  


                                                                     
source /public/home/mayuan/intel/oneapi/setvars.sh                       
ulimit -s unlimited                                                      
cd $PBS_O_WORKDIR                                                        
killall -9 pw.x                                                          


                                                                     
mpirun -np 28 /public/home/mayuan/software/qe-7.1/bin/pw.x -npool 4 <relax.in> relax.out 
check symmetry ops is consistent or not after vc-relax                   
grep "Sym. Ops." relax.out                                               
awk '/Begin final coordinates/,/End final coordinates/{print $0}' relax.out 
