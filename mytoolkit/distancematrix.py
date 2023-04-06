#!/usr/bin/env python

import sys
import numpy as np
print("-----------------------------------------HELP------------------------------------------")
print("You can use this scripts by example:")
print("distancematrix.py '2.8, 1.8, 1.1' ")
print("'2.8, 1.8, 1.1' represents the value of RCORE in POTCAR, its unit is Bohr unit !!!")
print("-----------------------------------------HELP------------------------------------------")
print("\n")
# bohr -> Angstrom
rcore   = np.array([eval(sys.argv[1])]) * 0.529 * 0.7 
dists_m = rcore + rcore.T
for row in dists_m:
    print(np.round(row[0],decimals=3), np.round(row[1],decimals=3), np.round(row[2],decimals=3),)
