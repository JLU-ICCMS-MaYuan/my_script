#!/bin/sh                                                                
#SBATCH  --job-name=ph_split                                             
#SBATCH  --output=log.ph_split.out.%j                                    
#SBATCH  --error=log.ph_split.err.%j                                     
#SBATCH  --partition=xieyu                                               
#SBATCH  --nodes=1                                                       
#SBATCH  --ntasks=48                                                     
#SBATCH  --ntasks-per-node=48                                            
#SBATCH  --cpus-per-task=1                                               


                                                                     
source /work/env/intel2018                                               
ulimit -s unlimited                                                      


                                                                     
mpirun -n 48 /work/software/q-e-qe-6.8/bin/pw.x -npool 4 <scf.fit.in> scf.fit.out
mpirun -n 48 /work/software/q-e-qe-6.8/bin/pw.x -npool 4 <scf.in> scf.out         
mpirun -n 48 /work/software/q-e-qe-6.8/bin/ph.x -npool 4 <split_ph.in> split_ph.out 
