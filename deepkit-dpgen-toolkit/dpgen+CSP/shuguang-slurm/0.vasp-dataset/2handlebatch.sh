#!/bin/bash

rm -fr {1..20}

for i in {1..20}; do 

mkdir $i
cp POSCAR_$i $i/POSCAR
cp POSCAR_$i $i/POSCAR_$i
cp INCAR_* $i/
cp POTCAR $i/
cp vasp.slurm $i/ 
done

wait
for i in {1..20}; do
cd $i
sbatch vasp.slurm
cd ..
done

