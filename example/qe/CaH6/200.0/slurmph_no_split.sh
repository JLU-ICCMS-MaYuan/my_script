#!/bin/sh                                                                                 
#SBATCH  --job-name=ph_no_split                                                           
#SBATCH  --output=log.ph_no_split.out                                                     
#SBATCH  --error=log.ph_no_split.err                                                      
#SBATCH  --partition=xieyu                                                                
#SBATCH  --nodes=1                                                                        
#SBATCH  --ntasks=48                                                                      
#SBATCH  --ntasks-per-node=48                                                             
#SBATCH  --cpus-per-task=1                                                                


                                                                                      
source /work/env/intel2018                                                                
ulimit -s unlimited                                                                       


                                                                                      
mpirun -n 48 /work/software/q-e-qe-6.8/bin/ph.x -npool 4 <ph_no_split.in> ph_no_split.out 
