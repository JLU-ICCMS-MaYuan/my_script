#!/usr/bin/env python3
import sys
import numpy 
from pymatgen.core.structure import Molecule
from pymatgen.symmetry.analyzer import PointGroupAnalyzer

filename = sys.argv[1]

mol = Molecule.from_file(filename)
sym_mol = PointGroupAnalyzer(mol).symmetrize_molecule()['sym_mol']


pg = PointGroupAnalyzer(sym_mol)

pgops = pg.get_symmetry_operations()
print(pg.get_pointgroup())

print(len(pgops))
for idx, op in enumerate(pgops):
    t = [str(int(vs)) for vs in op.translation_vector]
    r = [list(map(str, vs)) for vs in numpy.round(op.rotation_matrix, 2)]
    string = f'T{idx+1}=' + '{{{' + ', '.join(r[0]) + '}, ' + '{' + ', '.join(r[1]) + '}, ' + '{' + ', '.join(r[2]) + '}}, ' + '{' + ', '.join(t) + '}};'
    print(string)