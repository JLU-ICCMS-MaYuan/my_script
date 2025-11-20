'''
python supercell.py ./POSCAR [2,2,4]
'''
import sys
from pymatgen.core.structure import Structure
from pymatgen.io.vasp import Poscar

path = sys.argv[1]
supercell = eval(sys.argv[2])
struct = Structure.from_file(path)
struct.make_supercell(supercell)
Poscar(struct).write_file("SPOSCAR")