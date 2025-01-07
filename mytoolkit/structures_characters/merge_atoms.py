#!/usr/bin/env python3
import sys

from pymatgen.core.structure import Structure
from pymatgen.io.vasp import Poscar

filename = sys.argv[1]
tolerance = float(sys.argv[2])
struct = Structure.from_file(filename)

struct.merge_sites(tol=tolerance, mode='average')
Poscar(structure=struct).write_file("new.vasp")
