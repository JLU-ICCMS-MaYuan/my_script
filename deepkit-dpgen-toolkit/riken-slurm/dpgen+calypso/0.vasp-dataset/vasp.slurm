#!/bin/sh
#------ slurm option --------#
#SBATCH --partition=mpc_l
#SBATCH --account=hp240139
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --time=24:00:00

ulimit -s unlimited

#------- Program execution -------#
module load intelmpi/impi_23.2.0
module load intel/23.02.1


for i in 1 2 3 4 5; do
cp INCAR_$i INCAR
srun /lustre/home/h240012/soft/vasp.6.1.0/bin/vasp_std > vasp.log_$i 2>&1
cp CONTCAR POSCAR
done
