import os
import re
import shutil
from pathlib import Path
import logging
from argparse import ArgumentParser

from ase.io import read
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.vasp import Poscar

from qebin import qe_source_libs
from config import config

logger = logging.getLogger("qe_base")

class qe_base:
    def __init__(
        self,
        work_path: str,
        press: int,
        submit_job_system: str,
        input_file_path: Path,
    ) -> None:
        '''
        initial the qe_base, only different quantity the parameters
        work_path: decide 3 thing:
            1. the position directory filename/press
            2. the position of ppdirectory

        input_file_path:
            1. if .vasp: then the program will create directory filename/press,
               it will be named work_underpressure

            2. if relax.out: then the program will not create directory filename/press under the work_path,
               the parents of input_file_path will be work_underpressure

        if you do "relax":
            you must specify: 1. work_path, 2. press, 3. submit_job_system, 4. input_file_path
        if you do "scf":
            you must specify: 1. work_path,           3. submit_job_system  4. input_file_path
        '''
        self.work_path         = work_path
        self.press             = press          
        self.submit_job_system = submit_job_system
        self.input_file_path   = Path(input_file_path)
        
        if ".vasp" in self.input_file_path.name:
            self.input_file_name = self.input_file_path.name.split(".")[0]
            self.work_underpressure= Path(self.work_path).joinpath(self.input_file_name, str(self.press))
            if not self.work_underpressure.exists():
                self.work_underpressure.mkdir(parents=True, exist_ok=True)
        elif "relax.out" in self.input_file_path.name:
            self.work_underpressure= Path(self.input_file_path).parent
        else:
            self.work_underpressure= Path(self.work_path).joinpath(str(self.press))
            if not self.work_underpressure.exists():
                self.work_underpressure.mkdir(parents=True, exist_ok=True)

        self.ase_type          = read(self.input_file_path)
        self.struct_type       = AseAtomsAdaptor.get_structure(self.ase_type)
        self.get_struct_info(self.struct_type, self.work_underpressure)
        
        ############################ prepare pp directory #########################
        logger.info(f"create potcar dir in {self.work_path}")
        self.workpath_pppath = Path(self.work_path).joinpath("pp")
        if not self.workpath_pppath.exists():
            self.workpath_pppath.mkdir(parents=True)
        # 准备赝势 
        self.get_USPP(self.workpath_pppath)   
        ############################# done pp directory ##########################        
        
    def get_struct_info(self, struct, output_poscar):
        '''
        input parameter:
            struct:         pymatgen class structure
            output_poscar:  the output path of primitive cell
        return:
            None
        '''
        spa = SpacegroupAnalyzer(struct)
        # bstruct = spa.get_conventional_standard_structure()
        pstruct = spa.get_primitive_standard_structure()
        Poscar(pstruct).write_file(output_poscar.joinpath("PPOSCAR"))

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

    def get_USPP(self, workpath_pppath):
        '''
        input parameter:
            pseudopotention path 
        
        function(功能):
            if the target uspp doesn't exist:
                copy the target uspp that is in the potlib to custom potlib
            else:
                add the path of target uspp to the `self.final_choosed_pp`
        '''
        pp_files = os.listdir(workpath_pppath)
        self.final_choosed_pp = []
        for species in self.species:
            species_name = species.name.lower()
            dst_pps = filter(lambda file: re.search("^"+species_name+"\_", file.lower()), pp_files)
            dst_pps = list(dst_pps)
            if len(dst_pps) == 1:
                print(1)
                self.final_choosed_pp.append(dst_pps[0])
            elif len(dst_pps) == 0:
                logger.info(f"to prepare {species_name} uspp! ")
                dst_pp = self.get_single_uspp(species_name, workpath_pppath)
                self.final_choosed_pp.append(dst_pp)
            else:
                logger.error(f"find many POTCARs {dst_pps}")
            logger.info(f"species_name {species_name}")

    def get_single_uspp(self, species_name, workpath_pppath):
        '''
        function:
            obtain a single uspp 
        input:
            species_name: element name from pymatgen.composition.species.name
            workpath_pppath: 
        '''
        qe_USPP = os.path.abspath(qe_source_libs)
        if os.path.exists(qe_USPP):
            ppfiles = os.listdir(qe_USPP)
        else:
            raise FileExistsError("You may not set the potcar_source_libs, you can set it in `qebin.py` ")
        targetppfiles = filter(lambda file: re.search("^"+species_name+"\_", file.lower()), ppfiles)
        targetppnames = [pp for pp in targetppfiles]
        choosed_flag  = False
        while not choosed_flag:
            choosed_pp = input(f"{targetppnames}, \nplease input you want one\n")
            if choosed_pp in targetppnames:
                src_pp = Path(qe_USPP).joinpath(choosed_pp)
                dst_pp = Path(workpath_pppath).joinpath(choosed_pp)
                shutil.copy(src_pp, dst_pp)
                choosed_flag = True
            else:
                choosed_flag = False
            
        return choosed_pp