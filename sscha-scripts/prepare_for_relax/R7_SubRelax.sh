#!/bin/sh                           
#------ slurm option --------#
#SBATCH --partition=mpc
#SBATCH --account=hp240139
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=24:00:00

ulimit -s unlimited

#------- Program execution -------#
module load intelmpi/impi_23.2.0
module load intel/23.02.1

conda activate sscha
python Relax.py > OUT12.dat
