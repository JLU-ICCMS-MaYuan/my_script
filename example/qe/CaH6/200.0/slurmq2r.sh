#!/bin/sh                                                                  
#SBATCH  --job-name=q2r                                                    
#SBATCH  --output=log.q2r.out                                         
#SBATCH  --error=log.q2r.err                                          
#SBATCH  --partition=xieyu                                                 
#SBATCH  --nodes=1                                                         
#SBATCH  --ntasks=48                                                       
#SBATCH  --ntasks-per-node=48                                              
#SBATCH  --cpus-per-task=1                                                 


                                                                       
source /work/env/intel2018                                                 
ulimit -s unlimited                                                        


                                                                       
mpirun -n 48 /work/software/q-e-qe-6.8/bin/q2r.x -npool 4 <q2r.in> q2r.out 
grep nqs q2r.out > nqs                                                     
