import os
import re
import sys
import shutil
from pathlib import Path
from typing import List

import numpy as np
from ase.io import read
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.vasp import Poscar

from vasp.vaspbin import potcar_source_libs


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
        self.press             = int(press)          
        self.submit_job_system = submit_job_system
        self.input_file_path   = input_file_path
        self.mode              = mode

        self.input_file_name   = self.input_file_path.name.split(".")[0]
        self.ase_type          = read(self.input_file_path)
        try:
            self.struct_type       = AseAtomsAdaptor.get_structure(self.ase_type)
        except:
            print("\nNote: ------------------------------")
            print("    Read inputfile {} occurs something wrong!".format(self.input_file_name))
            sys.exit(1)
            
        print("Step 1 ------------------------------")
        print("    Create work_path under the work_path.")
        if self.work_path is None:
            self.work_path = Path.cwd()
            print(f"    You didn't specify the work_path, the default work_path is the current path {self.work_path}!")
        else:
            print(f"    You specify the work_path, the work_path will be '{self.work_path}' !")
            self.work_path = Path(self.work_path)
            if not self.work_path.exists():
                self.work_path.mkdir(parents=True)

        print("Step 2 ------------------------------")
        print("    Prepare the POSCAR and its coresponding primitive and conventional cell.")
        self.get_struct_info(self.struct_type, self.work_path)

        print("Step 3 ------------------------------")
        print("    Prepare the directory of `potcar_lib` and merge POTCARs ")
        print(f"   Pick up POTCARs from: \n        {potcar_source_libs}")
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
        # 获得倒格子
        self.reciprocal_plattice = self.get_reciprocal_lattice()

    def get_reciprocal_lattice(self):
        """
        获得一个胞的倒格子
        """
        a1 = self.cell_parameters[0]
        a2 = self.cell_parameters[1]
        a3 = self.cell_parameters[2]

        b1 = 2*np.pi*np.cross(a2, a3) / np.dot(a1, np.cross(a2, a3))
        b2 = 2*np.pi*np.cross(a3, a1) / np.dot(a2, np.cross(a3, a1))
        b3 = 2*np.pi*np.cross(a1, a2) / np.dot(a3, np.cross(a1, a2))
 
        self.reciprocal_plattice = np.array([b1, b2, b3])
        return self.reciprocal_plattice

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
                print(f"find many POTCARs {dst_pps}")
                sys.exit(1)
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
            choosed_pot_number = int(input("please input you want one(Enter the previous number, such as 1,2,3...)\n"))
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
            print("Merge POTCAR of element {}".format(elements))
            src_potcar_path_list = []
            for ele in elements:
                for pot in final_choosed_potcar:
                    if ele == pot.name:
                        src_potcar_path_list.append(pot)
            dst_potcar = dst_potcar_path.joinpath("POTCAR")
            print(f"Order is:")
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

    def write_evenly_kpoints(self, lattice, kspacing, kpoints_path):
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

        kpoints = Path(kpoints_path).joinpath("KPOINTS")
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
        print(f"lattice type: \n{ltype}")
        pstring = ltype.special_path
        print(f"the high symmetry points path: \n{pstring}")

        _plist  = [[ p for p in pp if not p.isdigit()] for pp in pstring.split(",")]
        print(f"the high symmetry points path: \n{_plist}")

        print(
            "please input the mode you want, just even input Number like 1 or 2\n",
            "0:  all_points:\n",
            "1:  first_group_points\n",
            "2:  second_group_points\n",
            "n:  n^Th_group_points\n",
            "1 2: first_group_points and second_group_points"
            "...."
            "Nothing to input, directly press ENTER, the default is all_points\n"
            )
        try:
            high_symmetry_type = list(map(int, input().split())) #将输入的整数字符串按照空格划分成列表并分别转化为整数类型并再转化为列表
        except:
            print("what you input is not an integer number, So use the `0:  all_points`")
            high_symmetry_type = [0]

        path_name_list = []
        if "," in pstring:
            if 0 in high_symmetry_type:
                path_name_list = list(chain.from_iterable(_plist))
                print(f"the choosed high symmetry points path is \n {path_name_list}")
            elif 0 not in high_symmetry_type:
                # path_name_list = max(_plist, key=len)
                for hst in high_symmetry_type:
                    path_name_list.extend(_plist[hst-1])
                print(f"the choosed high symmetry points path is \n {path_name_list}")
        else:
            path_name_list = [ pp for pp in pstring]
        
        # 这里的高对称点是指：这个倒空间的布里渊区有几种高对称点，并不是某一个路径的高对称点的列表。
        # 如果想获得某一个路径下高对称点的坐标，需要按照path_name_list的顺序依次获得相应的坐标。
        special_points   = ltype.get_special_points()
        path_coords      = [list(special_points[pname]) for pname in path_name_list]
        


        # 处理高对称点路径
        print("Print Fractional Coordinates of Reciprocal Lattice ! ")
        for name, coord in zip(path_name_list, path_coords):
            print("{}      {:<8.6f} {:<8.6f} {:<8.6f}".format(name, coord[0], coord[1], coord[2]))
            # < 表示左对齐，8.6f 表示留出8个位置并保留小数点后6位。



        print("The reciprocal lattice (without multiplating `unit_reciprocal_axis`)")
        for vector in self.reciprocal_plattice:
            print("{:<6.3f} {:<6.3f} {:<6.3f} ".format(vector[0], vector[1], vector[2]))



        print("Print projected high symmetry path")
        print("倒格子的单位是: 2pi/埃")
        path_name_coords = list(zip(path_name_list, path_coords))
        projected_path_name_coords = [[path_name_coords[0][0], path_name_coords[0][1][0]]]
        total_dist = 0
        for idx in range(1, len(path_name_coords)):
            current_name   = path_name_coords[idx][0]
            current_coords = np.dot(self.reciprocal_plattice, path_name_coords[idx][1])
            last_coords    = np.dot(self.reciprocal_plattice, path_name_coords[idx-1][1])
            dist = np.linalg.norm(current_coords-last_coords, 2)
            total_dist += dist
            projected_path_name_coords.append([current_name, total_dist])
        string_names = '  '.join(coord[0] for coord in projected_path_name_coords)
        string_coord = '  '.join(str(np.round(coord[1], 6)) for coord in projected_path_name_coords)
        print(string_names)
        print(string_coord)

        
        return path_name_list, path_coords

    def read_hspp(self, hsppfile_path):
        with open(hsppfile_path, "r") as f:
            lines = f.readlines()

        hspplist         = [line.strip() for line in lines[4:] if line != "\n"]
        path_name_coords = [hspplist[i] for i in range(0, len(hspplist), 2)] + [hspplist[-1]]
        print(path_name_coords)
        projected_path_name_coords = [[path_name_coords[0].split()[-1], list(map(float, path_name_coords[0].split()[:-1]))[0]]]
        total_dist = 0
        for idx in range(1, len(path_name_coords)):
            current_name   = path_name_coords[idx].split()[-1]
            current_coords = np.dot(self.reciprocal_plattice, list(map(float, path_name_coords[idx].split()[:-1])))
            last_coords    = np.dot(self.reciprocal_plattice, list(map(float, path_name_coords[idx-1].split()[:-1])))
            dist = np.linalg.norm(current_coords-last_coords, 2)
            total_dist += dist
            projected_path_name_coords.append([current_name, total_dist])
        string_names = '  '.join(coord[0] for coord in projected_path_name_coords)
        string_coord = '  '.join(str(np.round(coord[1], 6)) for coord in projected_path_name_coords)
        print(string_names)
        print(string_coord)

    def write_highsymmetry_kpoints(self, ase_type, kpoints_path):

        kpoints_filepath = kpoints_path.joinpath("KPOINTS")
        print("\nNote: --------------------------------")
        vaspkitflag = input("If you have installed vaspkit and you want to use it, input: Yes\n")
        if vaspkitflag:
            cwd = os.getcwd()
            os.chdir(kpoints_path)
            os.system('echo -e "3\n303" | vaspkit')
            shutil.copy("KPATH.in", "KPOINTS")
            os.system("sed -i '2s/.*/200/' KPOINTS")
            os.chdir(cwd)
        else: 
            path_name_list, path_coords = self.get_hspp(ase_type)
            pair_two_names  = [[path_name_list[i], path_name_list[i+1]] for i in range(len(path_name_list)-1)]
            pair_two_coords = [[path_coords[i], path_coords[i+1]] for i in range(len(path_coords)-1)]
            with open(kpoints_filepath, "w") as kp:
                kp.write("KPATH\n")
                kp.write("200\n")
                kp.write("Line-Mode\n")
                kp.write("Reciprocal\n")
                for two_names, two_coords in zip(pair_two_names, pair_two_coords):
                    kp.write("{:<10.8f} {:<10.8f} {:<10.8f} ! ${}$\n".format(two_coords[0][0],two_coords[0][1],two_coords[0][2], two_names[0]))
                    kp.write("{:<10.8f} {:<10.8f} {:<10.8f} ! ${}$\n".format(two_coords[1][0],two_coords[1][1],two_coords[1][2], two_names[1]))
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
        self.press             = int(press)          
        self.submit_job_system = submit_job_system
        self.input_file_path   = input_file_path
        self.mode              = mode
        
        self.input_file_name   = self.input_file_path.name.split(".")[0]
        if self.work_path is None:
            print("\nNote: -------------------")
            print("    You didn't specify the work_path, and this is parameters that must be set!!! the program will exit")
            sys.exit(1)
        else:
            self.work_path= Path(self.work_path).joinpath(str(self.press), self.input_file_name)
            if not self.work_path.exists():
                self.work_path.mkdir(parents=True)
            print("\nNote: --------------------")
            print("    Now {} will be created".format(self.work_path))

        self.ase_type          = read(self.input_file_path)
        self.struct_type       = AseAtomsAdaptor.get_structure(self.ase_type)
        self.get_struct_info(self.struct_type, self.work_path)
        
        ############################ prepare pp directory #########################
        self.workpath_pppath = Path(self.work_path).parent.parent.joinpath("potcar_lib")
        print(f"\nNote: --------------------")
        print(f"   To create potcar dir in \n {self.workpath_pppath}, it's convenient for you to choose POTCARs once.")
        if not self.workpath_pppath.exists():
            self.workpath_pppath.mkdir()
        # 准备赝势 
        self.get_potcar(self.work_path)   
