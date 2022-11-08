#!/bin/sh 

for a in `seq -w 16 18` 
do
cd $a
    sbatch 3step.sh
    cd ..
done
