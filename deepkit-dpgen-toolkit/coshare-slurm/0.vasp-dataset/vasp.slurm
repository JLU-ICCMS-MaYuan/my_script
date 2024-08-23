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

ulimit -s unlimited

for i in 1 2 3 4 5; do
cp INCAR_$i INCAR
srun /public/home/liuhanyu/software/apps/vasp.6.3.2/bin/vasp_std > vasp.log_$i 2>&1
cp CONTCAR POSCAR
done
