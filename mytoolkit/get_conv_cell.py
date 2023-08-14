#!/usr/bin/env python3

import sys
from pathlib import Path
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.vasp import Poscar

try:
    print("Please input filename and symprec, for example:")
    filename = sys.argv[1]
    symprec  = float(sys.argv[2])
except:
    print("You input nothing, so the default values will be used")
    print("    filename is POSCAR, symprec is 1e-3")
    filename = "POSCAR"

struct = Structure.from_file(filename)
spgaly = SpacegroupAnalyzer(struct, symprec = 0.001)
conv = spgaly.get_conventional_standard_structure()
prim = spgaly.get_primitive_standard_structure()

newfilename = Path(filename).parent.joinpath("conv-" + Path(filename).name)
Poscar(conv).write_file(newfilename)
