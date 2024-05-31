#!/usr/bin/env python3

import sys

from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

try:
    print("Please input filename and symprec, for example:")
    print("    get_spg.py CONTCAR 1e-5")
    filename = sys.argv[1]
    symprec  = float(sys.argv[2])
except:
    print("You input nothing, so the default values will be used")
    print("    filename is POSCAR, symprec is 1e-3")
    filename = "POSCAR"
    symprec  = 1e-3

struct = Structure.from_file(filename)
spgaly = SpacegroupAnalyzer(struct, symprec = 0.001)
symbol = spgaly.get_space_group_symbol()
spgnum = spgaly.get_space_group_number()

print("{} {}".format(symbol, spgnum))