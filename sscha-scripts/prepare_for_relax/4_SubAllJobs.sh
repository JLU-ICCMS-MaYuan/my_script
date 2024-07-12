#!/bin/bash
mkdir run_calculation
cd run_calculation
#mkdir run_calculation/tmp
mkdir tmp
mv *.pwo tmp
for i in `seq 1 100`
do
sbatch sub_$i.sh
sleep 3
done
cd ../

