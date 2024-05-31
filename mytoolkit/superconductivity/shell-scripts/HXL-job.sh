#!/bin/sh 

for a in `seq -w 4 8`
do
cd $a
    sbatch 3step.sh
    cd ..
done
