#!/bin/bash
#SBATCH  --job-name=opt_fine                                             
#SBATCH  --output=opt_fine.out.%j                                        
#SBATCH  --error=opt_fine.err.%j                                         
#SBATCH  --partition=lhy
#SBATCH  --nodes=1                                                       
#SBATCH  --ntasks=48                                                
#SBATCH  --ntasks-per-node=48                                            
#SBATCH  --cpus-per-task=1   

source /work/env/intel2018

ulimit -s unlimited
export I_MPI_ADJUST_REDUCE=3
export MPIR_CVAR_COLL_ALIAS_CHECK=0
conda activate /work/home/may/miniconda3/envs/deepmd

dp train input.json
