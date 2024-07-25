#!/bin/bash

numberoftasks=$1
echo "Total number of tasks " $numberoftasks
mkdir run_calculation
cd run_calculation
#mkdir run_calculation/tmp
mkdir tmp
mv *.pwo tmp

for i in `seq 1 $numberoftasks`
do
    sbatch sub_$i.sh
    sleep 3
done
cd ../
