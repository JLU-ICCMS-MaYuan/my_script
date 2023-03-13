from argparse import ArgumentParser
import os
import re
import logging
import shutil
from pathlib import Path
from typing import List

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
    create pressure directory (named work_underpressure) on the basis of work_path and press
    Create pseudopotential directory (named workpath_pppath) in the work_underpressure
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
        if self.work_path is None:
            self.work_underpressure = Path.cwd()
            self.work_path = self.work_underpressure.parent
        else:
            self.work_underpressure= Path(self.work_path).joinpath(str(int(self.press)))
            if not self.work_underpressure.exists():
                self.work_underpressure.mkdir(parents=True)

        self.ase_type          = read(self.input_file_path)
        self.struct_type       = AseAtomsAdaptor.get_structure(self.ase_type)
        self.get_struct_info(self.struct_type, self.work_underpressure)
        
        ############################ prepare pp directory #########################
        logger.info(f"create potcar dir in {self.work_underpressure}")
        self.workpath_pppath = Path(self.work_underpressure).joinpath("potcar_lib")
        if not self.workpath_pppath.exists():
            self.workpath_pppath.mkdir()
        # 准备赝势 
        self.get_potcar(self.work_underpressure)   
        ############################# done pp directory ##########################

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
        Poscar(pstruct).write_file(output_poscar.joinpath("PPOSCAR"))
        Poscar(bstruct).write_file(output_poscar.joinpath("BPOSCAR"))
        Poscar(struct).write_file(output_poscar.joinpath("POSCAR"))

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
                logger.info(f"to prepare {species_name} POTCAR! ")
                dst_pp = self.get_single_potcar(species_name)
                self.final_choosed_potcar.append(dst_pp)
            else:
                logger.error(f"find many POTCARs {dst_pps}")
        logger.info(f"the choosed pp is {self.final_choosed_potcar}")
        if len(self.final_choosed_potcar) == len(self.species):
            logger.info("merge")
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
            choosed_pot = input(f"{targetpotcarnames}, \nplease input you want one\n")
            if choosed_pot in targetpotcarnames:
                src_pp = Path(potcar_dir).joinpath(choosed_pot , "POTCAR")
                dst_pp = Path(self.workpath_pppath).joinpath(species_name)
                shutil.copy(src_pp, dst_pp)
                with open(Path(self.workpath_pppath).joinpath("potcar.log"), "a") as log:
                    log.write("{:<5}  {:<5}".format(species_name, choosed_pot))
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
        if self.work_underpressure.joinpath("POSCAR").exists():
            src_poscar = self.work_underpressure.joinpath("POSCAR") 
            elements = open(src_poscar, "r").readlines()[5].split()
            logger.info("The program will merge the POTCAR of element {}".format(elements))
            src_potcar_path_list = []
            for ele in elements:
                for pot in final_choosed_potcar:
                    if ele in pot.name:
                        if pot.exists():  
                            src_potcar_path_list.append(pot)
            dst_potcar = dst_potcar_path.joinpath("POTCAR")
            logger.info(f"POTCAR merge order is: {src_potcar_path_list}")
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
            self.work_underpressure = Path.cwd().joinpath(str(self.press), self.input_file_name)
            self.work_path = self.work_underpressure.parent
        else:
            self.work_underpressure= Path(self.work_path).joinpath(str(self.press), self.input_file_name)
            if not self.work_underpressure.exists():
                self.work_underpressure.mkdir(parents=True)

        self.ase_type          = read(self.input_file_path)
        self.struct_type       = AseAtomsAdaptor.get_structure(self.ase_type)
        self.get_struct_info(self.struct_type, self.work_underpressure)
        
        ############################ prepare pp directory #########################
        logger.info(f"create potcar dir in {self.work_path}")
        self.workpath_pppath = Path(self.work_path).joinpath("potcar_lib")
        if not self.workpath_pppath.exists():
            self.workpath_pppath.mkdir()
        # 准备赝势 
        self.get_potcar(self.work_underpressure)   
