#!/bin/bash
source /work/env/intel2020
a=$(grep -n 'Broadening   0.0200' gam.lines | cut -d ':' -f1)
b=$(grep -n 'Broadening   0.0250' gam.lines | cut -d ':' -f1)
sed -n ''$a','$b'p' gam.lines >>bb.dat
lastnum=`wc -l bb.dat | cut -d' ' -f1`
sed -i "$lastnum"d bb.dat
sed -i '1d' bb.dat
PREFIX=`grep prefix scf.in | cut -d "'" -f 2`
sed -i 's/lasih8.freq/'$PREFIX'.freq/g' gam.f90
ifort -o gam.x gam.f90
./gam.x
wait
cat bands.dat | awk '{print $1"  "$2}' >> "$PREFIX"_pho.dat

