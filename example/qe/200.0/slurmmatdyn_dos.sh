#!/bin/sh                                                                
##SBATCH  --job-name=matgen_dos                                           
##SBATCH  --output=log.matgen_dos.out                                     
##SBATCH  --error=log.matgen_dos.err                                      
##SBATCH  --partition=normal                                                  
##SBATCH  --nodes=1                                                       
##SBATCH  --ntasks=48                                                     
##SBATCH  --ntasks-per-node=48                                            
##SBATCH  --cpus-per-task=1                                               


                                                                     
source /work/home/may/intel/oneapi/setvars.sh --force                                                                
#ulimit -s unlimited                                                      


                                                                     
#mpirun -n 48 /work/home/may/software/qe-7.1/bin/matdyn.x -npool 4 <matdyn.dos.in> matdyn.dos.out          

mpirun -n 4 /work/home/may/software/qe-7.1/bin/matdyn.x <matdyn.dos.in> matdyn.dos.out          
