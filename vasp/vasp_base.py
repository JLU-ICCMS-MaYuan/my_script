import os
import re
import logging
import shutil
from pathlib import Path
from pprint import pprint

from ase.io import read
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.vasp import Poscar

from config import config

logger = logging.getLogger("vaspbase")

class vaspbase:

    def __init__(self, args) -> None:
        self._config = config(args).read_config()

        self.work_path         = self._config['work_path']
        self.press             = self._config["press"]          
        self.work_underpressure     = Path(self.work_path).joinpath(str(self.press))
        if not self.work_underpressure.exists():
            self.work_underpressure.mkdir()
        if "input_file_path" in self._config:
            self.input_file_path   = self._config["input_file_path"]
            self.ase_type          = read(self.input_file_path)
            self.struct_type       = AseAtomsAdaptor.get_structure(self.ase_type)
            self.get_struct_info(self.struct_type)
        else:
            logger.warning("you didn't specify the inputfile")
        self.work_path         = self._config["work_path"]      
        self.submit_job_system = self._config["submit_job_system"]
        
        if "workpath_pppath" in self._config: 
            self.workpath_pppath   = self._config["workpath_pppath"]
        else:
            self.workpath_pppath   = None


        if self.input_file_path is None:
            logger.error("please specify the inputfile *.vasp or POSCAR")
            raise ValueError ("please specify the inputfile *.vasp or POSCAR")
        ############################ prepare pp directory #########################
        logger.info(f"create potcar dir in {self.work_underpressure}")
        if self.workpath_pppath is None:
            self.workpath_pppath = Path(self.work_underpressure).joinpath("potcar_lib")
        if not self.workpath_pppath.exists():
            self.workpath_pppath.mkdir()
        # 准备赝势 
        self.get_potcar()   
        ############################# done pp directory ##########################

    def get_struct_info(self, struct):
        
        spa = SpacegroupAnalyzer(struct)
        # bstruct = spa.get_conventional_standard_structure()
        pstruct = spa.get_primitive_standard_structure()
        Poscar(pstruct).write_file(self.work_underpressure.joinpath("POSCAR"))

        # 处理PPOSCAR的pymatgen对象
        # 获得元素名称 和 每种元素的原子个数
        self.composition        = pstruct.composition.get_el_amt_dict()
        self.species            = pstruct.composition.elements
        # 获得体系的 化学配比
        self.system_name        = pstruct.composition.formula.replace(" ", "")
        # 获得元素种类的个数
        self.species_quantity   = len(self.composition)
        # 获得每种元素的相对原子质量
        self.all_atoms_quantity = int(sum(self.composition.values()))
        # 获得晶格矩阵
        self.cell_parameters    = pstruct.lattice.matrix
        # 获得原子分数坐标
        self.fractional_sites   = pstruct.sites

    def get_potcar(self):

        # 准备分离的赝势
        potlib_files = os.listdir(self.workpath_pppath)
        if len(potlib_files) != len(self.species):
            potcar_dir1 = os.path.abspath("/work/home/may/POT/vasp_pot1/potpaw_PBE54")
            potcar_dir2 = os.path.abspath("/public/home/mayuan/POT/vasp_pot1/potpaw_PBE54")
            if os.path.exists(potcar_dir1):
                potcar_dir = potcar_dir1
                potcarfiles = os.listdir(potcar_dir1)
            elif os.path.exists(potcar_dir2):
                potcar_dir = potcar_dir2
                potcarfiles = os.listdir(potcar_dir2)
            self.final_choosed_potcar = []
            for species in self.species:
                species_name = species.name
                logger.info(f"species_name {species_name}")
                targetpotcarfiles = filter(lambda file: re.search("^"+species_name, file), potcarfiles)
                targetpotcarnames = [pp for pp in targetpotcarfiles]
                choosed_flag  = False
                while not choosed_flag:
                    choosed_pot = input(f"{targetpotcarnames}, \nplease input you want one\n")
                    if choosed_pot in targetpotcarnames:
                        src_pp = os.path.join(potcar_dir, choosed_pot , "POTCAR")
                        dst_pp = os.path.join(self.workpath_pppath, "POTCAR-"+choosed_pot)
                        shutil.copy(src_pp, dst_pp)
                        choosed_flag = True
                        self.final_choosed_potcar.append("POTCAR-"+choosed_pot)
                    else:
                        choosed_flag = False
        else:
            self.final_choosed_potcar = potlib_files

        logger.info(f"the choosed pp is {self.final_choosed_potcar}")

        # 合并分立的POTCAR
        if not os.path.exists(self.workpath_pppath):
            raise FileExistsError("potcar_lib doesn't exist!")
        if self.work_underpressure.joinpath("POSCAR").exists():
            src_poscar = self.work_underpressure.joinpath("POSCAR") 
            elements = open(src_poscar, "r").readlines()[5].split()
            logger.info("element {}".format(elements))
            src_potcar_path_list = []
            for ele in elements:
                for pot in self.final_choosed_potcar:
                    if ele in pot:
                        src_potcar = self.workpath_pppath.joinpath(pot)
                        if src_potcar.exists():  
                            src_potcar_path_list.append(src_potcar)
            logger.info(f"{src_potcar_path_list}")
            dst_potcar = self.work_underpressure.joinpath("POTCAR")
            # 将多个POTCAR写入总的POTCAR中
            f = open(dst_potcar, "w")
            for potcar_indi in src_potcar_path_list:
                for line in open(potcar_indi):
                    f.writelines(line)
            f.close()
            os.system("grep TITEL {}".format(dst_potcar))
            print("完成了一个POTCAR\n\n\n")
        else:
            print("POSCAR not exist")




        
        
         

