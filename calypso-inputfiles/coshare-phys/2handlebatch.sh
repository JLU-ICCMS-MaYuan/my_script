#!/bin/bash

rm -fr {1..300}

for i in {1..300}; do 

mkdir $i
cp POSCAR_$i $i/POSCAR
cp INCAR_* $i/
cp POTCAR $i/
cp vasp.slurm $i/ 
done

for i in {1..300}; do
cd $i
sbatch vasp.slurm
cd ..
done

