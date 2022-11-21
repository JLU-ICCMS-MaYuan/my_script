#!/bin/sh                                     
#SBATCH  --job-name=calypso
#SBATCH  --output=log.calypso.out.%j             
#SBATCH  --error=log.calypso.err.%j              
#SBATCH  --partition=normal                       
#SBATCH  --nodes=1                            
#SBATCH  --ntasks=64
#SBATCH  --ntasks-per-node=64             
#SBATCH  --cpus-per-task=1                    

source /public/env/mpi_intelmpi-2021.3.0.sh
source /public/env/compiler_intel-compiler-2021.3.0.sh

ulimit -s unlimited
#export I_MPI_ADJUST_REDUCE=3
#export MPIR_CVAR_COLL_ALIAS_CHECK=0
export I_MPI_FABRICS=shm
export MKL_DEBUG_CPU_TYPE=5

rm -fr {1..300}

for i in {1..300}; do 

mkdir $i
cp POSCAR_$i $i/POSCAR
cp POSCAR_$i $i/POSCAR_$i
cp INCAR_* $i/
cp POTCAR $i/
# cp vasp.slurm $i/ 
done

for i in {1..300}; do
cd $i
for j in {1..4};do
cp INCAR_$j INCAR
timeout 14400s mpirun -n 64 /public/software/apps/vasp/intelmpi/5.4.4/bin/vasp_std > vasp.log 2>&1
cp CONTCAR POSCAR
done
cd ..
done

