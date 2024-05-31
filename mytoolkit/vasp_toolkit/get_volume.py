#!/usr/bin/env python3

import sys

from ase.io import read

try:
    filename = sys.argv[1]
except:
    filename = "POSCAR"
s = read(filename)
print(filename, s.get_volume())