#!/usr/bin/env python

import sys
import os
import re
import random 
from pathlib import Path

numberofstructures = eval(sys.argv[1])

with open("vsave", "w") as f:
    for i in range(1, numberofstructures+1):
        line = re.findall(r"\d+", os.popen(f"sed -n '6p' POSCAR_{i}").read())
        all_number = sum(list(map(int, line)))
        random_number = random.randint(1, 1000)
        f.write("{:<10} {:<10}\n".format(all_number, random_number))