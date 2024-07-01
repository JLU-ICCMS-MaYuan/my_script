#!/usr/bin/env python3

import sys

from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.vasp import Poscar

try:
    print("Please input filename and symprec, for example:")
    print("    get_spg.py CONTCAR 1e-5")
    filename = sys.argv[1]
    symprec  = float(sys.argv[2])
    struct = Structure.from_file(filename)
    spgaly = SpacegroupAnalyzer(struct, symprec = symprec)
    symbol = spgaly.get_space_group_symbol()
    spgnum = spgaly.get_space_group_number()
    conv_cell = spgaly.get_conventional_standard_structure()
    prim_cell = spgaly.get_primitive_standard_structure()
    Poscar(prim_cell).write_file("PPOSCAR")
    Poscar(conv_cell).write_file("BPOSCAR")
    print("{:<10} {:<10} {:<10}".format("symprec", "symbol", "spgnum"))
    print("{:<10.9f} {:<10} {:<10}".format(symprec, symbol, spgnum))
except:
    print("You input structures only, so the default values will be used")
    print("    filename is POSCAR, symprec is 1e-3")
    filename = sys.argv[1]
    struct = Structure.from_file(filename)
    print("{:<10} {:<10} {:<10}".format("symprec", "symbol", "spgnum"))
    for symprec in [0.1, 0.01, 0.001, 0.0001, 0.00001, 0.000001]:
        spgaly = SpacegroupAnalyzer(struct, symprec = symprec)
        symbol = spgaly.get_space_group_symbol()
        spgnum = spgaly.get_space_group_number()
        print("{:<10.6f} {:<10} {:<10}".format(symprec, symbol, spgnum))
