import os
import re
import shutil
import logging
import numpy as np
from pathlib import Path

from argparse import ArgumentParser

from ase.io import read
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.vasp import Poscar

from qe.qebin import qe_source_libs
from qe.config import config

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
            1. if relax.out: then the program will not create directory filename/press under the work_path,
               the parents of input_file_path will be work_underpressure
            2. if others:  then the program will create directory `/press`, it will be named work_underpressure
            ####3. if None  :  then input_file_path = Path(self.work_path).joinpath("relax.out")
        if you do "relax":
            you must specify: 1. work_path, 2. press, 3. submit_job_system, 4. input_file_path
        if you do "scf":
            you must specify: 1. work_path,           3. submit_job_system  4. input_file_path
        '''
        self.work_path         = work_path
        self.press             = press          
        self.submit_job_system = submit_job_system
        self.input_file_path   = Path(input_file_path)

        # 设置underpressure
        if ("relax.out" in self.input_file_path.name) or ("scf.out" in self.input_file_path.name) or ("scffit.out" in self.input_file_path.name):
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
        self.workpath_pppath = Path(self.work_underpressure).joinpath("pp")
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
        bstruct = spa.get_conventional_standard_structure()
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
        # 获得倒格矢
        # 注意，这是用于固体物理的标准倒数晶格，因数为2π
        self.reciprocal_plattice = self.get_reciprocal_lattice(self.work_underpressure.joinpath("scffit.out"))
        # self.reciprocal_plattice = pstruct.lattice.reciprocal_lattice
        # self.reciprocal_blattice = bstruct.lattice.reciprocal_lattice
        # 返回晶体倒数晶格，即没有2π的因子。
        # self.reciprocal_lattice_crystallographic = pstruct.lattice.reciprocal_lattice_crystallographic


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
            species_name = species.name
            dst_pps = get_pps_for_a_element(species_name, pp_files)
            if len(dst_pps) > 1:
                # 如果目标元素已经有多个赝势放在pp目录中. 例如S元素的三个赝势都放在了pp目录中，你需要手动选择哪个赝势用作计算
                choosed_flag  = False
                while not choosed_flag:
                    choosed_pp = input(f"The pp directory has the {dst_pps}! please input you want one\n")
                    if choosed_pp in dst_pps:
                        choosed_flag = True
                        self.final_choosed_pp.append(choosed_pp)
                    else:
                        choosed_flag = False
                print(self.final_choosed_pp)
            elif len(dst_pps) == 1:
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
            pp_files = os.listdir(qe_USPP)
        else:
            raise FileExistsError("You may not set the potcar_source_libs, you can set it in `qebin.py` ")

        targetppfiles    = get_pps_for_a_element(species_name, pp_files)
        targetppnames    = [pp for pp in targetppfiles]
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

    def get_reciprocal_lattice(self, scffit_out: Path):
        """
        读取scffit.out文件获得其中的reciprocal lattice
        Read the scffit.out file to get the reciprocal lattice
        """

        if not scffit_out.exists():
            logger.warning("You have run scffit for getting `scffit.out`")
            return None
            
        reciprocal_lattice = []
        for i in range(1,4):
            b = os.popen(f"sed -n '/b({i})/p' {scffit_out}").read()
            b = re.findall(r"-?\d+\.\d+", b)
            b  = list(map(float, b))
            reciprocal_lattice.append(b)
        return np.array(reciprocal_lattice)


def get_pps_for_a_element(species_name:str, pp_files:list):
    """
    species_name: 元素名
    pp_files    : 存储着赝势名的列表
    该函数可以从存储着赝势名的列表中挑选出指定元素的全部赝势
    支持识别的赝势名称格式分别为：
        H.pbe-van_bm.UPF.txt  首字母大写的元素名.****
        h_pbe_v1.uspp.F.UPF   首字母小写的元素名_****
    """
    # 这里写的实在是太精妙了。我用2个花括号套在一起解决了正则表达式传入花括号的问题
    namelenth = len(species_name)
    if namelenth == 1:
        qepp_name_patter = re.compile(f'''^[{species_name.capitalize()}|{species_name.lower()}]{'{1,1}'}[\.|\_]''')
    elif namelenth == 2:
        qepp_name_patter = re.compile(f'''^[{species_name.capitalize()}|{species_name.lower()}]{'{2,2}'}[\.|\_]''')
    else:
        logger.warning(f"Inconceivable that the length of the element name is not 1 or 2, it is {namelenth}, its name is {species_name}")
    dst_pps          = list(filter(lambda file: qepp_name_patter.findall(file), pp_files))

    return dst_pps