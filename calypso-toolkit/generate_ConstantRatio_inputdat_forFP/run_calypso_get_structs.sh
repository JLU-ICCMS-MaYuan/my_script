#!/bin/bash

for i in `cat fixed_comp.name`; do 
    cd $i
    ./calypso.x
    cd ../
done
