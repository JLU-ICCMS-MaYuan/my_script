#!/usr/bin/env python3

from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

mesh_input = input("请输入网格点(例如:10 10 10)千万不要在数字间加入逗号")
mesh_input = [ int(mesh) for mesh in mesh_input.split() ]
shift_input= input("请输入平移量(例如:0.5 0.5 0.5)千万不要在数字间加入逗号")
shift_input= [ float(shift) for shift in shift_input.split() ]
stc = Structure.from_file("POSCAR")
spg = SpacegroupAnalyzer(stc, symprec=0.01)

ir_reduced_kpoints = spg.get_ir_reciprocal_mesh(mesh=mesh_input, is_shift=shift_input)
for coord, weight in ir_reduced_kpoints:
    print("{:<8.6f} {:<8.6f} {:<8.6f} {:<4}".format(coord[0], coord[1], coord[2], weight))

print("不可约k点个数是{}".format(len(ir_reduced_kpoints)))

