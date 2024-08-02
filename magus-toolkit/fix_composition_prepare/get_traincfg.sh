#!/bin/bash

for i in {00..36}; do 
    cd $i/labeling  
    for j in fp*; do
        cd $j
        mlp convert-cfg OUTCAR train.cfg --input-format=vasp-outcar 
        cat train.cfg >> ../../../train.cfg
        cd ..
    cd ../..
 done
