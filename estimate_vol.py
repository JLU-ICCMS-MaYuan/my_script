#!/usr/bin/env python

import sys
import numpy as np
print("-----------------------------------------HELP------------------------------------------")
print("You can use this scripts by example:")
print("estimate.py '2.8, 1.8, 1.1' '2, 4, 20'")
print("'2.8, 1.8, 1.1' represents the value of RCORE in POTCAR, its unit is Bohr unit !!!")
print("'2, 4, 20' respresents the value of NUMBER of every specie")
print("-----------------------------------------HELP------------------------------------------")
print("\n")

rcore = np.array(eval(sys.argv[1]))
num   = np.array(eval(sys.argv[2]))

if len(rcore) != len(num):
    raise ValueError("rcore != num")

v = [4*np.power(r*0.529, 3)*np.pi*n/3 for r, n in zip(rcore, num)]
volume = np.asarray(v).sum(axis=0)
print("volume=%s"%(volume))

