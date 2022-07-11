#!/work/home/may/miniconda3/bin/python3
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
                opt_dir = os.path.join(substitute_sourcedir, pstruct_copycopy.composition.formula.replace(" ", ""))
                if not os.path.exists(opt_dir):
                    os.makedirs(opt_dir)
                Poscar(pstruct_copycopy).write_file(os.path.join(opt_dir, name))
                Poscar(pstruct_copycopy).write_file(os.path.join(opt_dir, "POSCAR"))


    # 准备分立的POTCAR
    all_element = M_element + X_element + ["H"]
    POT_dir = os.path.abspath("/work/home/may/POT/vasp_pot1/potpaw_PBE54")
    if pot_individual and subtitute_structure_directory is not None:
        # 提前准备好一个赝势库的空目录
        pot_lib = os.path.join(subtitute_structure_directory, "potcar_lib")
        if not os.path.exists(pot_lib):
            os.makedirs(pot_lib)
        for element in all_element:
            print("\n\n ------------element={}-------------".format(element))
            dirs_list = os.listdir(POT_dir)
            for d in dirs_list:
                if element in d:
                    print(d)
            input_valid_flag = False
            while not input_valid_flag:
                potcar_dir = input("请选择你要使用的POTCAR: \n")
                if potcar_dir in dirs_list:
                    src_potcar =os.path.join(POT_dir, potcar_dir, "POTCAR")
                    dst_potcar = os.path.join(pot_lib, element)
                    if os.path.exists(src_potcar):
                        shutil.copy(src_potcar, dst_potcar)
                        worklog = os.path.join(pot_lib, "pot_selected.log")
                        with open(worklog, "a") as w:
                            w.write("{} {}\n".format(potcar_dir, src_potcar))
                        input_valid_flag = True
                    else:
                        input_valid_flag = False
                else:
                    input_valid_flag = False

    # 合并分立的POTCAR
    if pot_merge and subtitute_structure_directory is not None:
        pot_lib = os.path.join(subtitute_structure_directory, "potcar_lib")
        pot_name_list = os.listdir(pot_lib)
        if not os.path.exists(pot_lib):
            raise FileExistsError("potcar_lib doesn't exist!")
        for root, dirs, files in os.walk(subtitute_structure_directory):
            for f in files:
                if re.search(r"\.vasp$", f):
                    src_poscar = os.path.join(root, f)
                    elements = open(src_poscar, "r").readlines()[5].split()
                    input("element {}".format(elements))
                    src_potcar_path_list = []
                    for ele in elements:
                        if ele in pot_name_list:
                            src_potcar = os.path.join(pot_lib, ele)
                            if os.path.exists(src_potcar):  
                                src_potcar_path_list.append(src_potcar)
                    input("src_potcar_path_list={}".format(src_potcar_path_list))
                    dst_potcar = os.path.join(root, "POTCAR")
                    f = open(dst_potcar, "w")
                    for potcar_indi in src_potcar_path_list:
                        for line in open(potcar_indi):
                            f.writelines(line)
                    f.close()

                    os.system("grep TITEL {}".format(dst_potcar))
                    input("完成了一个POTCAR")
                    

                    

                    

                    



 
