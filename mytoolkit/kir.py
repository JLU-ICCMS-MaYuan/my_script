#!/usr/bin/env python
import math

from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer


def get_irreduced_k(mesh_input):

    shift_input= [0, 0, 0]
    stc = Structure.from_file("POSCAR")
    spg = SpacegroupAnalyzer(stc, symprec=0.01)

    ir_reduced_kpoints = spg.get_ir_reciprocal_mesh(mesh=mesh_input, is_shift=shift_input)
    # for coord, weight in ir_reduced_kpoints:
    #     print("{:<8.6f} {:<8.6f} {:<8.6f} {:<4}".format(coord[0], coord[1], coord[2], weight))

    print("总的k点是{}".format(mesh_input[0]*mesh_input[1]*mesh_input[2]))
    print("不可约k点个数是{}".format(len(ir_reduced_kpoints)))

if __name__ == "__main__":
    n_ks = input("请输入k点的密度(例如: 5 5 5)\n").split()
    if len(n_ks) == 3:
        n_ks = [int(nk_i) for nk_i in n_ks]
        get_irreduced_k(n_ks)
    else:
        print("你输入的k点个数必须是三个, 且用空格分割")
