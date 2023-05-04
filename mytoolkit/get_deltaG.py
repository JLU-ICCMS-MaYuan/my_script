#!/usr/bin/env python

import sys
import math

print("Note: --------------------")
print("This script just can calculate deltaG with considering the entropy of configeration for Binary alloys")
dH = input("Please input enthalpy per atoms, if you don't know it, just press ENTER\n")
if not dH:
    dH = float(input("Please input enthalpy for a cell obtained from OUTCAR\n"))
    N  = float(input("Please input total number of atoms for a cell\n"))
    dH = dH / N
    print("So enthalpy per atom is {:<12.8f} eV/atom".format(dH))
else:
    dH = float(dH)
    print("So enthalpy per atom is {:<12.8f} eV/atom".format(dH))
xA = float(input("Please input percentage of A-element\n"))
xB = 1 - xA
print("So percentage of B-element is {:.3f}".format(xB))
kB = 8.6173324e-5 #eV/K

Sconf_PerAtom = - kB * (xA * math.log(xA) + xB * math.log(xB))
print("So entropy is {:<12.8f} eV/K".format(Sconf_PerAtom))
print("{:<10} {:<12}".format("T", "dG"))
for T in range(0, 5100, 100):
    dG = dH-T*Sconf_PerAtom 
    print("{:<10} {:<12.8f}".format(T, dG))
