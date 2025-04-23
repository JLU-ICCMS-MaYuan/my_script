import os
import re
import sys
import shutil
import logging
from pathlib import Path
from collections import Counter

import numpy as np

from ase.io import read
from ase.atoms import Atoms
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.vasp import Poscar

logger = logging.getLogger(__name__)

class epw_base:
    def __init__(
        self,
        work_path: str,
        submit_job_system: str,
        input_file_path: Path,
    ) -> None:
        '''
        epw_main.py -i 输入文件路径 -w 工作目录 -p 压强 -j 运行方式
        说明: 输入文件路径 和工作路径说明 : 
        1. 如果输入文件是relax.out, 那么工作目录会被强行限制为输入文件是relax.out的上一级路径。
        2. 如果输入文件是其它路径:
            1. 如果指定了-w,   那么就是在-w指定的路径下开展计算
            2. 如果没有指定-w, 那么就默认所有的计算都在当前指定epw_main.py命令的目录下运行. 并不会额外创建一个压强值命名的目录作为最终工作目录

        说明: 压强
        1. 指定结构优化的压强, 单位是GPa.

        说明: 运行方式
        1. bash 代表本地使用bash shell运行。
        2. slurm 代表使用slurm脚本运行。
        3. pbs 代表使用pbs脚本运行。

        说明: 赝势文件最终存放在压强命名的目录的下面
        '''
        self.work_path         = work_path 
        self.submit_job_system = submit_job_system
        self.input_file_path   = Path(input_file_path)

        # if ("relax.out" in self.input_file_path.name) or ("scf.out" in self.input_file_path.name) or ("scffit.out" in self.input_file_path.name):
            # self.work_path = Path(self.input_file_path).parent
        if self.work_path is None:
            self.work_path = Path.cwd()
        else:
            self.work_path = Path(self.work_path)
            if not self.work_path.exists():
                self.work_path.mkdir(parents=True, exist_ok=True)
        logger.info(f"work_path: {self.work_path.absolute()}")

        try:
            self.ase_type = read(self.input_file_path)
        except:
            logger.error("When reading `{}` file, the program get something wrong, you need to check it !!!".format(self.input_file_path))
            sys.exit(1)

        self.struct_type = AseAtomsAdaptor.get_structure(self.ase_type)
        self.get_struct_info(self.struct_type, self.work_path)
        ############################ prepare pp directory #########################
        logger.info(f"Create pp in {self.work_path.absolute()}")
        logger.info(f"Pick up UPF from: {epw_source_libs}")
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
        bstruct = spa.get_conventional_standard_structure()
        pstruct = spa.get_primitive_standard_structure()
        # Poscar(pstruct).write_file(output_poscar.joinpath("PPOSCAR"))
        # 备份原始结构
        logger.debug(f"finish back up origin inputed structure file into workpath: {self.work_path.absolute()}")
        if not self.work_path.joinpath(self.input_file_path.name).exists():
            # shutil.copy(self.input_file_path, self.work_path)
            backup_inputfile = self.work_path.joinpath(self.input_file_path.name)
            os.system(f"cp -f {self.input_file_path} {backup_inputfile}")

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
        self.cell_parameters    = self.get_cell(struct) 
        # 获得原子分数坐标
        self.fractional_sites   = self.get_coords(struct)
        # 获得倒格矢
        # 注意，这是用于固体物理的标准倒数晶格，因数为2π
        self.reciprocal_plattice= self.get_reciprocal_lattice(struct)
        # self.reciprocal_plattice = pstruct.lattice.reciprocal_lattice
        # self.reciprocal_blattice = bstruct.lattice.reciprocal_lattice
        # 返回晶体倒数晶格，即没有2π的因子。
        # self.reciprocal_lattice_crystallographic = struct.lattice.reciprocal_lattice_crystallographic
        self.celldm1 = self.get_celldm1()
