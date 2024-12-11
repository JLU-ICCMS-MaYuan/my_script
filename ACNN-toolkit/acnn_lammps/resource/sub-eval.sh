#!/bin/sh
#===========================================================
# Configure DSUB resources
#===========================================================
#DSUB --job_type cosched
#DSUB -n MMM
#DSUB -A root.jildxwlxyljlstdui
#DSUB -q root.default
#DSUB -R cpu=16
#DSUB -N 1
##DSUB -pn cccs-share-agent-[117,126,129]
#DSUB -oo donau.out.%J
#DSUB -eo donau.err.%J
#===========================================================
#/home/share/jildxwlxyljlstdui/home/lijx/playground/torchdemo-dev/build/acnn -train ./inv2 > train-log 2>&1
export PATH=/home/share/jildxwlxyljlstdui/home/lijx/software/torchdemo/build2:$PATH

module purge
export OMP_NUM_THREADS=16

module use /home/HPCBase/workspace/public/software/modules
module load compilers/gcc/gcc9.3.0
module load compilers/bisheng/bisheng2.1.0
module load mpi/hmpi/1.1.1/hmpi1.1.1-bisheng2.1.0

ldd $(which acnn)


emodel model-restart/model-10000  train_dt> ss 2>&1