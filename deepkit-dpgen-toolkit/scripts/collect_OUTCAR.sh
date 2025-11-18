#!/bin/bash

directories=("1.Lu-Li-H.[1-4][1-4][10-40]" "2.Lu-Li" "3.Li-H" "4.Lu-H" "5.H" "6.Li" "7.Lu")
i=0
mkdir group
for d in "${directories[@]}"; do
    for n in {1..500..1}; do
        i=$((i+1))
        cd ${d}/${n} && cp OUTCAR ../../group/OUTCAR_${i} && cd ../..
    done 
done
