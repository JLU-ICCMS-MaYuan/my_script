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

for i in {1..4}; do 
cp INCAR_$i INCAR
killall -9 vasp_std
sleep 3
timeout 14400s mpirun -np 48 /work/software/vasp.6.1.0/vasp_std > vasp.log_${i} 2>&1
wait
cp CONTCAR POSCAR
done

rm -rf CHG* WAVECAR

