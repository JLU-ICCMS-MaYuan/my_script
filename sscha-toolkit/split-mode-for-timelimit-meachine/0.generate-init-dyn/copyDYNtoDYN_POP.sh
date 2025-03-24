#!/bin/bash

total_dyns=$1
prefix=$2

if [ ! -d "1.relax" ]; then
    mkdir -p ../1.relax
fi 

for i in $(seq 1 $total_dyns); do
    cp ${prefix}.dyn$i ../1.relax/dyn_pop0_$i
done
