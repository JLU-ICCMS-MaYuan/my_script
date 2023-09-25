#!/bin/bash
# how to use:
# 	3handle.sh step_1 
#	that means to create a directory named step_1
mkdir $1
cp POSCAR_* $1/
cp CONTCAR_* $1/
cp OUTCAR_* $1/
