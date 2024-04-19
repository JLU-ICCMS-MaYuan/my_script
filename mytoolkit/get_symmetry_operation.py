import sys

from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.symmetry.analyzer import SpacegroupOperations
filename = sys.argv[1]

struct = Structure.from_file(filename)
spg = SpacegroupAnalyzer(struct)
pgops = spg.get_point_group_operations()
symOps = spg.get_symmetry_operations()
print(pgops)
print(len(pgops))
print(symOps)

