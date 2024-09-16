#!/usr/bin/env python3
import sys
from pymatgen.core.structure import Structure
from pymatgen.io.vasp import Poscar

inputfilename = sys.argv[1]
deleteelements = sys.argv[2].split()

struct = Structure.from_file(inputfilename)
struct.remove_species(deleteelements)

try:
    outputfilename =  sys.argv[3]
    Poscar(structure=struct).write_file(outputfilename)
except:
    Poscar(structure=struct).write_file("deleted_"+inputfilename)
