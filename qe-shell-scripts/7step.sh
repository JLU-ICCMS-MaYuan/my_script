#!/bin/sh
#SBATCH  --job-name=pwscf
#SBATCH  --output=log.out.%j
#SBATCH  --error=log.err.%j
#SBATCH  --partition=lhy
#SBATCH  --nodes=1
#SBATCH  --ntasks=48
#SBATCH  --ntasks-per-node=48
#SBATCH  --cpus-per-task=1

source /work/env/intel2018

srun hostname | sort | uniq >> /tmp/nodefile.$$
NP=`srun hostname | wc -l`
#mpirun -np 48 /work/home/hxl/opt/soft/qe/q-e-qe-6.8.new/q-e-qe-6.8/bin/lambda.x -npool 4 <lambda.in> lambda.out
mpirun -np 48 /work/home/hxl/opt/soft/qe/q-e-qe-6.8-600nmodex/bin/lambda.x -npool 4 <lambda.in> lambda.out
