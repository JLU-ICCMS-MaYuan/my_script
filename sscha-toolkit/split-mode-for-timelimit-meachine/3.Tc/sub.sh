#!/bin/bash

nq=`sed -n 2p Nb4H14.dyn0`
for i in `seq 2 $nq`; do
    cd $i
    rm slurm-*
    sbatch s5_PhAssignQ.sh
    cd ..
done