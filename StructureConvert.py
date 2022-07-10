#!/work/home/mayuan/miniconda3/envs/cage/bin/python3
'''
StructureFormatConvert.py -if  ICSD.cif -cif2poscar
     -if  ICSD.cif 
     -cif2poscar
'''
import os
import shutil
from argparse import ArgumentParser

from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.cif import CifWriter
from pymatgen.io.vasp import Poscar


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument(
        "-if",
        "--input-file",
        default=None,
        type=str,
        required=True,
        dest="input_file",
        help="请输入cif文件的目录"
    )
    parser.add_argument(
        "-cif2poscar",
        "--cif-convert-to-poscar",
        default=False,
        action="store_true",
        dest="cif_convert_to_poscar",
        help="是否将cif转化为poscar"
    )
    parser.add_argument(
        "-ma",
        "--merge-atoms",
        default=None,
        action="store",
        type=float,
        dest="merge_atoms",
        help="是否将cif转化为poscar, 精度几何"
    )
    args = parser.parse_args()

    input_file            = args.input_file
    cif_convert_to_poscar = args.cif_convert_to_poscar
    merge_atoms           = args.merge_atoms

    if input_file is not None:
        input_file = os.path.abspath(input_file); print(input_file); input()
        struct = Structure.from_file(input_file)
        if merge_atoms is not None:
            struct.merge_sites(merge_atoms, mode='average')
            Poscar(struct).write_file("MPOSCAR")
        bstruct = SpacegroupAnalyzer(struct).get_conventional_standard_structure()
        pstruct = SpacegroupAnalyzer(struct).get_primitive_standard_structure()
    
    if cif_convert_to_poscar and struct is not None:
        Poscar(pstruct).write_file("PPOSCAR")
        Poscar(bstruct).write_file("BPOSCAR")
        shutil.copy("PPOSCAR", "POSCAR")
    






