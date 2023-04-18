from argparse import ArgumentParser
import os
import re
import logging
import shutil
from pathlib import Path
from typing import List

import numpy as np
from ase.io import read
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.vasp import Poscar

from vasp.vaspbin import potcar_source_libs
from vasp.config import config

logger = logging.getLogger("vaspbase")

class vasp_base:
    """
    Read the basic information from ARGS. The basic information includes:
        work_path
        press
        input_file_path
        submit_job_system
    create pressure directory (named work_path) on the basis of work_path and press
    Create pseudopotential directory (named workpath_pppath) in the work_path
    """
    def __init__(
        self,
        work_path: str,
        press: int,
        submit_job_system: str,
        input_file_path: Path,
        mode: str,
    ) -> None:

        self.work_path         = work_path
        self.press             = press          
        self.submit_job_system = submit_job_system
        self.input_file_path   = input_file_path
        self.mode              = mode

        self.input_file_name   = self.input_file_path.name.split(".")[0]
        self.ase_type          = read(self.input_file_path)
        self.struct_type       = AseAtomsAdaptor.get_structure(self.ase_type)
        
        print("Step 1 ------------------------------")
        print("    Create work_path under the work_path.")
        print("    If you didn't specify the work_path, the default work_path is the current path and the work_path is its parent path !")
        print("    If you specify the work_path, the work_path will be 'work_path+(press/scf/eband/eledos)' !")
        if self.work_path is None:
            self.work_path = Path.cwd()
        else:
            self.work_path = Path(self.work_path)
            if not self.work_path.exists():
                self.work_path.mkdir(parents=True)

        print("Step 2 ------------------------------")
        print("    Prepare the POSCAR and its coresponding primitive and conventional cell.")
        self.get_struct_info(self.struct_type, self.work_path)

        print("Step 3 ------------------------------")
        print("    Prepare the directory of `potcar_lib` and merge POTCARs ")
        self.workpath_pppath = Path(self.work_path).joinpath("potcar_lib")
        if not self.workpath_pppath.exists():
            self.workpath_pppath.mkdir() 
        self.get_potcar(self.work_path)
        print(f"potcar_lib is in \n {self.workpath_pppath}")


    def get_struct_info(self, struct, output_poscar):
        """
        function: get information of the `struct`
        input parameter:
            struct: the class of the `pymatgen.core.structure.Structure`
            output_poscar: the output path of your structure
        """ 
        spa = SpacegroupAnalyzer(struct)
        # bstruct = spa.get_conventional_standard_structure()
        pstruct = spa.get_primitive_standard_structure()
        bstruct = spa.get_conventional_standard_structure()
        Poscar(pstruct).write_file(output_poscar.joinpath("cell-primitive.vasp"))
        Poscar(bstruct).write_file(output_poscar.joinpath("cell-conventional.vasp"))
        Poscar(struct).write_file(output_poscar.joinpath("POSCAR"))
        print("NOTES: ------------------------------ ")
        print("    You really confirm the inputfile, such as POSCAR, ***.vasp,  is what you want !")
        # 处理PPOSCAR的pymatgen对象
        # 获得元素名称 和 每种元素的原子个数
        self.composition        = struct.composition.get_el_amt_dict()
        self.species            = struct.composition.elements
        # 获得体系的 化学配比
        self.system_name        = struct.composition.formula.replace(" ", "")
        # 获得元素种类的个数
        self.species_quantity   = len(self.composition)
        # 获得每种元素的相对原子质量
        self.all_atoms_quantity = int(sum(self.composition.values()))
        # 获得晶格矩阵
        self.cell_parameters    = struct.lattice.matrix
        # 获得原子分数坐标
        self.fractional_sites   = struct.sites

    def get_potcar(self, dst_potcar_path: Path):
        """
        prepare your every single potcar, and then merge all potcars to a complete potcar
        input parameter:
            dst_potcar_path: it determines the directory position of the final complete potcar,
                            please note it's not the self.workpath_pppath
        """
        # 准备分离的赝势
        self.final_choosed_potcar = []
        for species in self.species:
            species_name = species.name
            potlib_files = os.listdir(self.workpath_pppath)
            patter = re.compile(r"^{}$".format(species_name))
            dst_pps = list(filter(lambda file: re.search(patter, file), potlib_files))
            if len(dst_pps) == 1:
                dst_pp = Path(self.workpath_pppath).joinpath(dst_pps[0])
                self.final_choosed_potcar.append(dst_pp)
            elif len(dst_pps) == 0:
                print(f"to prepare {species_name} POTCAR! ")
                dst_pp = self.get_single_potcar(species_name)
                self.final_choosed_potcar.append(dst_pp)
            else:
                logger.error(f"find many POTCARs {dst_pps}")
        if len(self.final_choosed_potcar) == len(self.species):
            self.merge_potcar(dst_potcar_path, self.final_choosed_potcar)

    def get_single_potcar(self, species_name):
        """
        according to the given name of element, copy the POTCAR to the potcar_lib
        input:
            species_name
        rerutn:
            the class Path of POTCAR of the given element.
        """
        potcar_dir = os.path.abspath(potcar_source_libs)
        if os.path.exists(potcar_dir):
            potcarfiles = os.listdir(potcar_dir)
        else:
            raise FileExistsError("You may not set the potcar_source_libs, you can set it in `vaspbin.py` ")
        targetpotcarfiles = filter(lambda file: re.search("^"+species_name, file), potcarfiles)
        targetpotcarnames = [pp for pp in targetpotcarfiles]
        choosed_flag  = False
        while not choosed_flag:
            for i, p in enumerate(targetpotcarnames):
                print(f"{i+1}.  ", p)
            choosed_pot_number = input("please input you want one(Enter the previous number, such as 1,2,3...)\n")
            choosed_pot = targetpotcarnames[choosed_pot_number-1]
            if choosed_pot in targetpotcarnames:
                src_pp = Path(potcar_dir).joinpath(choosed_pot , "POTCAR")
                dst_pp = Path(self.workpath_pppath).joinpath(species_name)
                shutil.copy(src_pp, dst_pp)
                with open(Path(self.workpath_pppath).joinpath("potcar.log"), "a") as log:
                    log.write("{:<10}  {:<10}\n".format(species_name, choosed_pot))
                choosed_flag = True
            else:
                choosed_flag = False
        return dst_pp

    def merge_potcar(
        self, 
        dst_potcar_path: Path,
        final_choosed_potcar: List[Path],
        ):
        """
        dst_potcar_path: specify the path of  merged POTCAR
        final_choosed_potcar: a list that saves all potcars that will be merged
        """
        # 合并分立的POTCAR
        # 第一步：判断potcar_lib是否存在
        # 第二步：找到目标POSCAR，并读取其元素组成放在 src_potcar_path_list
        # 第三步：将完成的POTCAR合并出来，并且合并到目标POSCAR所在的目录
        if not os.path.exists(self.workpath_pppath):
            raise FileExistsError("potcar_lib doesn't exist!")
        if self.work_path.joinpath("POSCAR").exists():
            src_poscar = self.work_path.joinpath("POSCAR") 
            elements = open(src_poscar, "r").readlines()[5].split()
            print("The program will merge the POTCAR of element {}".format(elements))
            src_potcar_path_list = []
            for ele in elements:
                for pot in final_choosed_potcar:
                    if ele in pot.name:
                        if pot.exists():  
                            src_potcar_path_list.append(pot)
            dst_potcar = dst_potcar_path.joinpath("POTCAR")
            print(f"POTCAR merge order is:")
            for path in src_potcar_path_list:
                print(path.name)
            # 将多个POTCAR写入总的POTCAR中 方法1
            # f = open(dst_potcar, "w")
            # for potcar_indi in src_potcar_path_list:
            #     for line in open(potcar_indi):
            #         f.writelines(line)
            # f.close()
            # os.system("grep TITEL {}".format(dst_potcar))
            # print("完成了一个POTCAR\n\n\n")
            # 将多个POTCAR写入总的POTCAR中 方法2
            os.system(f"rm -fr {str(dst_potcar)}")
            for potcar_indi in src_potcar_path_list:
                os.system(f"cat {str(potcar_indi)} >> {str(dst_potcar)}")
        else:
            print("POSCAR not exist")

    def write_evenly_kpoints(self, kspacing, position):
        lattice = self.cell_parameters
        a1 = lattice[0]
        a2 = lattice[1]
        a3 = lattice[2]
        b1 = np.cross(a2, a3) / np.dot(a1, np.cross(a2, a3))
        b2 = np.cross(a3, a1) / np.dot(a2, np.cross(a3, a1))
        b3 = np.cross(a1, a2) / np.dot(a3, np.cross(a1, a2))
        b1_mol = np.linalg.norm(b1)
        b2_mol = np.linalg.norm(b2)
        b3_mol = np.linalg.norm(b3)
        N1 = np.ceil(b1_mol*2*np.pi / float(kspacing))
        N2 = np.ceil(b2_mol*2*np.pi / float(kspacing))
        N3 = np.ceil(b3_mol*2*np.pi / float(kspacing))

        kpoints = Path(position).joinpath("KPOINTS")
        with open(kpoints, "w") as kp:
            kp.write("from kspacing is {}\n".format(kspacing))
            kp.write("0\n")
            kp.write("Gamma\n")
            kp.write("{:<5} {:<5} {:<5}\n".format(str(N1), str(N2), str(N3)))
            kp.write("{:<5} {:<5} {:<5}\n".format("0", "0", "0"))

    def create_kpoints_by_pymatgen(self, pmg_struct, output_kpoints, kdensity):
        """
        automatic_density_by_length方法 根据输入的结构(第一个参数), 各个维度k点的密度(第二个参数), 是否强制使用Gamma方法产生k点(第三个参数)
        重点说说如何设置各个方向的密度。该方法依据的公式是: 
            number of k points along this direction = density of k along this derection / parameter of this direction
        所以用户需要根据超胞的大小, 合理估计一个k点密度, 然后得到该方向的k点数
        """
        from pymatgen.io.vasp import Kpoints
        print("creat KPOINTS by `automatic_density_by_length`")
        kpoints = Kpoints.automatic_density_by_lengths(
            structure=pmg_struct, 
            length_densities=kdensity,
            force_gamma=True)
        kpoints.write_file(output_kpoints)
        print(kpoints)

    def get_hspp(self, ase_type):

        from itertools import chain

        ltype   = ase_type.cell.get_bravais_lattice()
        pstring = ltype.special_path
        _plist  = [[ p for p in pp if not p.isdigit()] for pp in pstring.split(",")]

        print(f"the high symmetry points path is {_plist}")

        print(
            "please input the mode you want, just even input Number like 1 or 2\n",
            "'1  all_points'\n",
            "'2  main_points'\n",
            "Nothing to input, directly press ENTER, the default is main_points\n"
            )
        high_symmetry_type = input()

        if not high_symmetry_type:
            high_symmetry_type = "2" # default
        if "," in pstring:
            if high_symmetry_type == "1":
                path_name_list = list(chain.from_iterable(_plist))
                print(f"the choosed high symmetry points path is \n {path_name_list}")
            elif high_symmetry_type == "2":
                path_name_list = max(_plist, key=len)
                print(f"the choosed high symmetry points path is \n {path_name_list}")
        else:
            path_name_list = [ pp for pp in pstring]
        
        # 这里的高对称点是指：这个倒空间的布里渊区有几种高对称点，并不是某一个路径的高对称点的列表。
        # 如果想获得某一个路径下高对称点的坐标，需要按照path_name_list的顺序依次获得相应的坐标。
        special_points   = ltype.get_special_points()
        path_coords      = [list(special_points[pname]) for pname in path_name_list]
        
        print(f"the choosed high symmetry path is\n")
        for name, coord in zip(path_name_list, path_coords):
            print("{}      {:<8.6f} {:<8.6f} {:<8.6f}".format(name, coord[0], coord[1], coord[2]))
            # < 表示左对齐，8.6f 表示留出8个位置并保留小数点后6位。
        return path_name_list, path_coords

    def write_highsymmetry_kpoints(self, ase_type, kpoints_path):
        path_name_list, path_coords = self.get_hspp(ase_type)
        pair_two_names  = [[path_name_list[i], path_name_list[i+1]] for i in range(len(path_name_list)-1)]
        pair_two_coords = [[path_coords[i], path_coords[i+1]] for i in range(len(path_coords)-1)]
        with open(kpoints_path, "w") as kp:
            kp.write("KPATH\n")
            kp.write("50\n")
            kp.write("Line-Mode\n")
            kp.write("Reciprocal\n")
            for two_names, two_coords in zip(pair_two_names, pair_two_coords):
                kp.write("{:<10.8f} {:<10.8f} {:<10.8f} ! ${}$\n".format(two_coords[0][0],two_coords[0][1],two_coords[0][2], two_names[0]))
                kp.write("{:<10.8f} {:<10.8f} {:<10.8f} ! ${}$\n".format(two_coords[1][0],two_coords[1][1],two_coords[1][2], two_names[0]))
                kp.write("\n")

class vaspbatch_base(vasp_base):

    def __init__(
        self,
        work_path: str,
        press: int,
        submit_job_system: str,
        input_file_path: Path,
        mode: str,
    ) -> None:

        self.work_path         = work_path
        self.press             = press          
        self.submit_job_system = submit_job_system
        self.input_file_path   = input_file_path
        self.mode              = mode
        
        self.input_file_name   = self.input_file_path.name.split(".")[0]
        if self.work_path is None:
            self.work_path = Path.cwd().joinpath(str(self.press), self.input_file_name)
            self.work_path = self.work_path.parent
        else:
            self.work_path= Path(self.work_path).joinpath(str(self.press), self.input_file_name)
            if not self.work_path.exists():
                self.work_path.mkdir(parents=True)

        self.ase_type          = read(self.input_file_path)
        self.struct_type       = AseAtomsAdaptor.get_structure(self.ase_type)
        self.get_struct_info(self.struct_type, self.work_path)
        
        ############################ prepare pp directory #########################
        print(f"create potcar dir in \n {self.work_path}")
        self.workpath_pppath = Path(self.work_path).joinpath("potcar_lib")
        if not self.workpath_pppath.exists():
            self.workpath_pppath.mkdir()
        # 准备赝势 
        self.get_potcar(self.work_path)   
