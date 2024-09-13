import sys

from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

filename = sys.argv[1]

s = Structure.from_file(filename)
print(s)
spg = SpacegroupAnalyzer(s)
spgops = spg.get_space_group_operations()
print(spg.get_space_group_symbol())

print(len(spgops))
for idx, op in enumerate(spgops):
    t = [str(int(vs)) for vs in op.translation_vector]
    r = [list(map(str, map(int, vs))) for vs in op.rotation_matrix]
    string = f'T{idx+1}=' + '{{{' + ', '.join(r[0]) + '}, ' + '{' + ', '.join(r[1]) + '}, ' + '{' + ', '.join(r[2]) + '}}, ' + '{' + ', '.join(t) + '}};'
    print(string)

