#!/usr/bin/env python

import sys
import numpy as np
rcore = np.array([eval(sys.argv[1]), eval(sys.argv[2]), eval(sys.argv[3])])
# bohr -> Angstrom
rcore = rcore * 0.529

DistanceOfIon = np.zeros((3,3))
for id, i in enumerate(rcore):
    row = i + rcore
    DistanceOfIon[id,:] = row

for row in DistanceOfIon:
    print(np.round(row[0],decimals=3), np.round(row[1],decimals=3), np.round(row[2],decimals=3),)
