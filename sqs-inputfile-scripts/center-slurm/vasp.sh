#!/bin/sh 
#SBATCH  --job-name=sqsvasp
#SBATCH  --output=log.out.%j
#SBATCH  --error=log.err.%j
#SBATCH  --partition=lhy
#SBATCH  --nodes=1
#SBATCH  --ntasks=48
#SBATCH  --ntasks-per-node=48
#SBATCH  --cpus-per-task=1


source /work/env/intel2018
#source /work/home/may/intel/oneapi/setvars.sh --force
export I_MPI_ADJUST_REDUCE=3
export MPIR_CVAR_COLL_ALIAS_CHECK=0
ulimit -s unlimited

srun hostname | sort | uniq >> /tmp/nodefile.$$
NP=`srun hostname | wc -l`

for i in {1..5}; do
    killall -9 vasp_std
    mpirun -np 48 /work/software/vasp.5.4.4/vasp_std_5.4.4 > vasp.log_$i 2>&1
    cp OUTCAR OUTCAR_$i
    cp CONTCAR POSCAR
    cp CONTCAR CONTCAR_$i
done
