#!/bin/bash
for i in `cat fixed_comp.name`; do 
    cd $i
    ./batchsubmit.sh
    cd ../
done
