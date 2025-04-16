#!/bin/bash

for i in 0 2 3 6 7; do 
    echo POSCAR-$i
    cd POSCAR-$i
    cp ../10convertoML.py .
    python 10convertoML.py OUTCAR POSCAR_$i.xsf
    mv step*.xsf ../datasets_La1_vs_Y3
    cd ..
done
