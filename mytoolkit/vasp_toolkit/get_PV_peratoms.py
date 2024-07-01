#!/usr/bin/env python3

import os
import sys
import re

try:
    PV = os.popen('grep "P V=" OUTCAR ' + "| awk -F= '{print $3}' | tail -n 1").read().strip('\n')
    if not PV:
        PV = 100000000000000.0
except:
    print("   OUTCAR有问题读不出来PV项")
    PV = 100000000000000.0

print("{:<12.8f}".format(float(PV)))