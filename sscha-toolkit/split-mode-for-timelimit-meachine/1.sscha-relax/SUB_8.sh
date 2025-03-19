#!/bin/bash
#SBATCH -J qe
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -p xahcnormal

module purge
module load compiler/intel/2021.3.0
module load mpi/intelmpi/2021.3.0
conda activate sscha
. 8_ReStart.sh
