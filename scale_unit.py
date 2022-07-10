#!/work/home/mayuan/miniconda3/envs/cage/bin/python3
import os
import numpy as np
from argparse import ArgumentParser
from pymatgen.core.structure import Molecule
from pymatgen.io.xyz import XYZ


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument(
        '-xyz',
        '--xyz-file',
        dest="xyz_file",
        default=None,
        action='store',
        type=str,
        help="指定xyz文件的位置"
    )
    parser.add_argument(
        '-c',
        '--centered-molecule',
        dest="centered_molecule",
        default=False,
        action='store_true',
        help="是否中心对称分子"
    )
    args = parser.parse_args()

    xyz_file          = args.xyz_file
    centered_molecule = args.centered_molecule

    xyz_file_path = os.path.abspath(xyz_file)

    if os.path.exists(xyz_file_path):

        mol     = Molecule.from_file(xyz_file_path)
        mol     = mol.get_centered_molecule()
        mol_box = mol.get_boxed_structure(a=30, b=30, c=30, no_cross=True)
        scaling_flag = True
        while scaling_flag:
            print(np.sort(mol_box.distance_matrix, axis=1)[...,1:5])
            scaling_factor = input("请输入放缩因子:") # 如果input什么也不输入，那么就是一个空的字符串
            if scaling_factor:
                old_volume = mol_box.volume
                new_volume = float(old_volume) * float(scaling_factor)
                mol_box.scale_lattice(new_volume)
            else:
                scaling_flag = False
       
        if centered_molecule:
            XYZ(mol_box).write_file(xyz_file_path)
            mol = Molecule.from_file(xyz_file_path)
            mol = mol.get_centered_molecule()
            XYZ(mol).write_file(xyz_file_path)
        else:
            XYZ(mol_box).write_file(xyz_file_path)







            