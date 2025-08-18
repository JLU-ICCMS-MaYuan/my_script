#!/bin/bash
#irrep -Ecut=50 -code=vasp -kpnames=GM -EF=15.0734  > out
#irrep -Ecut=100 -code=vasp -kpoints=200 -EF=15.0734  > out
irrep -Ecut=900 -code=vasp -kpoints=200 -kpnames=GM -EF=15.0734  > out
#irrep -Ecut=50 -code=vasp -kpoints=6 -kpnames=GM,H,N,P,GM,N -EF=13.4499  > out
#irrep -Ecut=50 -code=vasp -kpoints=1 -kpnames=GM,H,N,P -EF=13.4499 -IBstart=5 -IBend=10 > out
cp out myout
