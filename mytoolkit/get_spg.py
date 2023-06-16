#!/usr/bin/env python3

import sys

from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

filename = sys.argv[1]

struct = Structure.from_file(filename)
spgaly = SpacegroupAnalyzer(struct, symprec = 0.001,)
symbol = spgaly.get_space_group_symbol()
spgnum = spgaly.get_space_group_number()

print("{} {}".format(symbol, spgnum))