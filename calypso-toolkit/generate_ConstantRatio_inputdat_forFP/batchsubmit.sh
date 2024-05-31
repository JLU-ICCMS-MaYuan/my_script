#!/bin/bash

POSCAR_num=$(find . -maxdepth 1 -type f -name "POSCAR_*" | wc -l)

for ((i=1; i<=POSCAR_num; i++)); do 

mkdir $i
cp POSCAR_$i $i/POSCAR
cp POSCAR_$i $i/POSCAR_$i
cp INCAR_* $i/
cp POTCAR $i/
cp vasp.pbs $i/ 
done

for ((i=1; i<=POSCAR_num; i++)); do 
cd $i
qsub vasp.pbs
cd ..
done

