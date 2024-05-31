#!/usr/bin/env python3

import os
import sys
import re

try:
    dH = os.popen('grep "free  energy   TOTEN" OUTCAR' + " | tail -n 1 | awk '{print $ 5}'").read().strip('\n')
    if not dH:
        dH = 100000000000000.0
except:
    print("   OUTCAR有问题读不出来焓值")
    dH = 100000000000000.0

print("{:<12.8f}".format(float(dH)))