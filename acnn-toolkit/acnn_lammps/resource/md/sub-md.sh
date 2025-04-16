#!/bin/sh
#===========================================================
# Configure DSUB resources
#===========================================================
#DSUB --job_type cosched
#DSUB -n md-al
#DSUB -A root.jildxwlxyljlstdui
#DSUB -q root.default
#DSUB -R cpu=120
#DSUB -N 1
#DSUB -oo donau.out.%J
#DSUB -eo donau.err.%J
#===========================================================
#/home/share/jildxwlxyljlstdui/home/lijx/playground/torchdemo-dev/build/acnn -train ./inv2 > train-log 2>&1

module purge

module use /home/HPCBase/workspace/public/software/modules
module load compilers/gcc/gcc9.3.0
module load compilers/bisheng/bisheng2.1.0
module load mpi/hmpi/1.1.1/hmpi1.1.1-bisheng2.1.0

mpirun -np 120 ~/playground/torchdemo-dev/bin/lmp_mpi  -i in.lmp > 3-md-log 2>&1
#echo "OOS" >> 3-md-log
