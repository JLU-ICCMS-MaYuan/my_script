#!/bin/bash
source /work/env/intel2020
ifort -o dir2cart.x dir2cart.f90
ifort -o kpath_20.x kpath_20.f

./dir2cart.x > kpath.in
mv kpath.in kpath.ini
./kpath_20.x  
