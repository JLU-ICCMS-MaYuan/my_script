#!/bin/sh
#===========================================================
# Applying for resources
#===========================================================
#DSUB --job_type cosched
#DSUB -n vasp-exg
#DSUB -A root.jildxwlxyljlstdui
#DSUB -q root.default
#DSUB -R cpu=120
#DSUB -N 1
#DSUB -oo donau.out.%J
#DSUB -eo donau.err.%J

## -R  no. cores that a single host can apply for. Max: cpu=128;gpu=4
## -N  task copys (Nodes?)

#===========================================================
# load envs
#===========================================================
module purge
module use /home/HPCBase/workspace/public/software/modules
module load compilers/bisheng/bisheng2.5.0
module load libs/fftw/3.3.10/bisheng2.5.0-hmpi1.3.10
module load libs/openblas/openblas0.3.6_bisheng2.5.0
module load libs/scalapack/2.1.0/bisheng2.5.0_hmpi1.3.0
module load mpi/hmpi/1.3.1/hmpi1.3.1-bisheng2.5.0

export PATH=/home/HPCBase/workspace/public/software/apps/vasp.6.3.2/bin:$PATH
export PATH=/home/jildxwlxyljlstdui/cccs-share02/vasp/vasp.6.3.2-20240305/bin:$PATH

#mpirun -np 120 -x PATH=$PATH -x OMP_NUM_THREADS=1 vasp_std > log 2>&1
mpirun -np 120 -x PATH=$PATH -x OMP_NUM_THREADS=1 vasp_gam > log 2>&1




