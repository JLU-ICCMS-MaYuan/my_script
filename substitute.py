#!/work/home/mayuan/miniconda3/envs/cage/bin/python3
'''
检查POSCAR和POTCAR顺序是否匹配
根据POSCAR 合成 POTCAR
批量替换 POSCAR中某一个元素

使用例子：
    substitute.py -f LaBH10.vasp
'''

import os
import re
import shutil
from argparse import ArgumentParser
from itertools import combinations
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.vasp import Poscar
from pymatgen.core.periodic_table import Element



if __name__ == "__main__":

    parser = ArgumentParser()
    parser.add_argument(
        "-f",
        "--substituted-vasp-file",
        default=None,
        type=str,
        dest="substituted_vasp_file",
        help="请输入结构所在的路径"
    )
    parser.add_argument(
        "-d",
        "--subtitute-structure-directory",
        action="store",
        dest="subtitute_structure_directory",
        default=str,
        help="指定那些替换好的结构所在的母目录"
    )

    args = parser.parse_args()
    substituted_vasp_file         = args.substituted_vasp_file   
    subtitute_structure_directory = args.subtitute_structure_directory          # 准备POSCAR 建立运行opt的子目录
  
    M_element = [
        "Na", "K" , "Ru", "Cs",
        "Ca", "Sr", "Ba",
        "Sc", "Y" ,
        "La", "Ce", "Pr", "Tb",
        "Ac", "Th", "Pa",
        "Zr", "Hf"
        ]
    X_element = [
        "Li", "Be", "B", "C", "N",
        "Al", "Si", "P", "S"
    ]

    M_element_M_element = [
        
    ]

    # 准备POSCAR 和 POTCAR 建立运行opt的子目录
    if substituted_vasp_file is not None:
        src_filepath = os.path.abspath(substituted_vasp_file) # 未被移动到opt_dir中的vasp结构文件路径
        struct       = Structure.from_file(src_filepath)
        spa          = SpacegroupAnalyzer(struct)
        pstruct      = spa.get_primitive_standard_structure()

        substitute_sourcedir = os.path.abspath("substitute")
        if not os.path.exists(substitute_sourcedir):
            os.makedirs(substitute_sourcedir)
        # 替换金属元素
        for M_ele in M_element:
            pstruct_copy = pstruct.copy()
            pstruct_copy.replace_species({Element("La"):Element(M_ele)})
            # 替换轻元素
            for X_ele in X_element:
                pstruct_copycopy = pstruct_copy.copy()
                pstruct_copycopy.replace_species({Element("B"):Element(X_ele)})
                name = pstruct_copycopy.composition.formula.replace(" ", "")+".vasp"
                Poscar(pstruct_copycopy).write_file(os.path.join(substitute_sourcedir, name))
