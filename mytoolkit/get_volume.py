#!/usr/bin/env python3

import sys

from pymatgen.core.structure import Structure

try:
    filename = sys.argv[1]
except:
    filename = "POSCAR"
s = Structure.from_file(filename)
print(filename, s.volume)