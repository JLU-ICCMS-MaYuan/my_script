#!/bin/bash

for i in {1..300}; do
cp $i/POSCAR_$i POSCAR_$i
cp $i/CONTCAR CONTCAR_$i
cp $i/OUTCAR OUTCAR_$i
done

