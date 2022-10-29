#!/bin/sh                                                                
##SBATCH  --job-name=matgen                                               
##SBATCH  --output=log.matgen.out                                         
##SBATCH  --error=log.matgen.err                                          
##SBATCH  --partition=None                                                  
##SBATCH  --nodes=1                                                       
##SBATCH  --ntasks=48                                                     
##SBATCH  --ntasks-per-node=48                                            
##SBATCH  --cpus-per-task=1                                               

                                                                     
source /work/home/may/intel/oneapi/setvars.sh --force                                                                
ulimit -s unlimited                                                      

#mpirun -np 48 /work/home/may/software/qe-7.1/bin/matdyn.x -npool 4 <matdyn.in> matdyn.out                  
mpirun -np 1 /work/home/may/software/qe-7.1/bin/matdyn.x <matdyn.in_freq> matdyn.out_freq
