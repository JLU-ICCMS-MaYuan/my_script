#!/bin/sh                                                                
#SBATCH  --job-name=relax                                             
#SBATCH  --output=log.relax.out.%j                                        
#SBATCH  --error=log.relax.err.%j                                         
#SBATCH  --partition=xieyu                                               
#SBATCH  --nodes=1                                                       
#SBATCH  --ntasks=48                                                     
#SBATCH  --ntasks-per-node=48                                            
#SBATCH  --cpus-per-task=1                                               


                                                                     
source /work/env/intel2018                                               
ulimit -s unlimited                                                      


                                                                     
mpirun -n 48 /work/software/q-e-qe-6.8/bin/pw.x -npool 4 <relax.in> relax.out 
check symmetry ops is consistent or not after vc-relax                   
grep "Sym. Ops." relax.out                                               
awk '/Begin final coordinates/,/End final coordinates/{print $0}' relax.out 
