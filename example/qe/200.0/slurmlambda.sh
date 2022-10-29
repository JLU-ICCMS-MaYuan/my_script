#!/bin/sh                                                                
##SBATCH  --job-name=lambda                                           
##SBATCH  --output=log.lambda.out                                    
##SBATCH  --error=log.lambda.err                                     
##SBATCH  --partition=xieyu                                               
##SBATCH  --nodes=1                                                       
##SBATCH  --ntasks=48                                                     
##SBATCH  --ntasks-per-node=48                                            
##SBATCH  --cpus-per-task=1                                               


                                                                     
source /work/home/may/intel/oneapi/setvars.sh --force                                                                
ulimit -s unlimited                                                      

mpirun -np 1  /work/home/may/software/qe-7.1/bin/lambda.x <lambda.in> lambda.out 
