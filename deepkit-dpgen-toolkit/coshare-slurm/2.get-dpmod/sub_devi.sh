#!/bin/sh                           
#SBATCH  --job-name=mayqe                      
#SBATCH  --output=log.out                       
#SBATCH  --error=log.err                       
#SBATCH  --partition=intel6430
#SBATCH  --nodes=1                          
#SBATCH  --ntasks=64
#SBATCH  --ntasks-per-node=64
#SBATCH  --cpus-per-task=1                         

source /data/home/mayuan/intel/oneapi/setvars.sh

source /public/home/liuhanyu/workplace/mayuan/software/deepmd-kit/bin/activate 


python calypso_run_model_devi.py --all_models ../gen_stru_analy.000/graph.000.pb ../gen_stru_analy.000/graph.001.pb ../gen_stru_analy.000/graph.002.pb ../gen_stru_analy.000/graph.003.pb --type_map Ce Sc H
