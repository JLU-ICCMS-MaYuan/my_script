#!/bin/sh                                                                
#SBATCH  --job-name=scf                                                  
#SBATCH  --output=log.scf.out.%j                                             
#SBATCH  --error=log.scf.err.%j                                              
#SBATCH  --partition=xieyu                                               
#SBATCH  --nodes=1                                                       
#SBATCH  --ntasks=48                                                     
#SBATCH  --ntasks-per-node=48                                            
#SBATCH  --cpus-per-task=1                                               


                                                                     
source /work/env/intel2018                                               
ulimit -s unlimited                                                      


                                                                     
mpirun -n 48 /work/software/q-e-qe-6.8/bin/pw.x -npool 4 <scf.in> scf.out
