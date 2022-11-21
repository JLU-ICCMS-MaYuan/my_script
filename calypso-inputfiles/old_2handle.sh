#!/bin/bash

for i in {1..300}; do
cp $i/CONTCAR CONTCAR_$i
cp $i/OUTCAR OUTCAR_$i
done

