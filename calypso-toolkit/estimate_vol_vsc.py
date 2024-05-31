#!/usr/bin/env python

import sys
import numpy as np
print("-----------------------------------------HELP------------------------------------------")
print("-----------------------------Estimate volume for VSC-----------------------------------")
print("You can use this scripts by example: estimate_vol_vsc.py 'r1, r2, r3' ")
print("    eg: estimate.py '2.8, 1.8, 1.1' ")
print("'2.8, 1.8, 1.1' represents the value of RCORE in POTCAR, its unit is Bohr unit !!!")
print("-----------------------------------------HELP------------------------------------------")
print("\n")

rcore = np.array(eval(sys.argv[1]))
v = [4*np.power(r*0.529, 3)*np.pi/3.0 for r in rcore]
volume = 1.3*np.asarray(v).sum(axis=0)/len(rcore)
print("volume/atom=%.3f"%(volume))

