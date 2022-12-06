#!/bin/sh                                     


rm -fr {1..300}

for i in {1..300}; do 

mkdir $i
cp POSCAR_$i $i/POSCAR
cp POSCAR_$i $i/POSCAR_$i
cp INCAR_* $i/
cp POTCAR $i/
# cp vasp.slurm $i/ 
done

