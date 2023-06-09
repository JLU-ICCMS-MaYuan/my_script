#!/usr/bin/env python3

import os
import sys
import re

try:
    dH = os.popen("grep enthalpy OUTCAR | tail -n 1 | awk '{print $ 5}'").read().strip('\n')
except:
    print("   OUTCAR有问题读不出来焓值")
    dH = 100000000000000.0

try:
    begin_id = os.popen('grep -n "position of ions in cartesian coordinates" OUTCAR').read().split(":")[0]
    N = 0; row_id=int(begin_id)
    while True:
        row_id = row_id+1
        content  = os.popen("sed -n '{}p' OUTCAR".format(row_id)).read().strip('\n')
        corrds   = re.findall(r"[-+]?\d+\.\d+", content)
        if len(corrds) == 3:
            N += 1
        else:
            break
except:
    N = 1

print("{:<12.8f} {:<3}".format(float(dH)/N, N))