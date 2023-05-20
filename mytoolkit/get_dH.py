#!/usr/bin/env python3

import os
import sys

try:
    dH = os.popen("grep enthalpy OUTCAR | tail -n 1 | awk '{print $ 5}'").read().strip('\n')
except:
    print("   OUTCAR有问题读不出来焓值")
    sys.exit(1)

begin_id = os.popen('grep -n "position of ions in cartesian coordinates" OUTCAR').read().split(":")[0]
N = 0; row_id=int(begin_id)
while True:
    row_id = row_id+1
    content  = os.popen("sed -n '{}p' OUTCAR".format(row_id)).read().strip().split()
    if len(content) == 3:
        N += 1
    else:
        break

print("{:<12.8f} {:<3}".format(float(dH)/N, N))
