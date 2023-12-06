import os
import re
import sys
import time
import shutil
import logging
import numpy as np
from pathlib import Path


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
        qe_main.py -i 输入文件路径 -w 工作目录 -p 压强 -j 运行方式
        说明: 输入文件路径 和工作路径说明 : 
        1. 如果输入文件是relax.out, 那么工作目录会被强行限制为输入文件是relax.out的上一级路径。
        2. 如果输入文件是其它路径:
            1. 如果指定了-w,   那么就是在-w指定的路径下开展计算
            2. 如果没有指定-w, 那么就默认所有的计算都在当前指定qe_main.py命令的目录下运行. 并不会额外创建一个压强值命名的目录作为最终工作目录

        说明: 压强
        1. 指定结构优化的压强, 单位是GPa.

        说明: 运行方式
        1. bash 代表本地使用bash shell运行。
        2. slurm 代表使用slurm脚本运行。
        3. pbs 代表使用pbs脚本运行。

        说明: 赝势文件最终存放在压强命名的目录的下面
        '''
        self.work_path         = work_path
        self.press             = press          
        self.submit_job_system = submit_job_system
        self.input_file_path   = Path(input_file_path)

        if ("relax.out" in self.input_file_path.name) or ("scf.out" in self.input_file_path.name) or ("scffit.out" in self.input_file_path.name):
            self.work_path = Path(self.input_file_path).parent
        elif self.work_path is None:
            self.work_path = Path.cwd()
        else:
            self.work_path = Path(self.work_path)
            if not self.work_path.exists():
                self.work_path.mkdir(parents=True, exist_ok=True)
                
        try:
            self.ase_type = read(self.input_file_path)
        except:
            print("\nNote: --------------------")
            print("    When reading `{}` file, the program get something wrong, you need to check it !!!".format(self.input_file_path))
            sys.exit(1)

        self.struct_type = AseAtomsAdaptor.get_structure(self.ase_type)
        self.get_struct_info(self.struct_type, self.work_path)
        ############################ prepare pp directory #########################
        print("\nNote: ----------")
        print(f"    Create pp in {self.work_path}")
        print(f"    Pick up UPF from: \n        {qe_source_libs}")
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
        Poscar(pstruct).write_file(output_poscar.joinpath("PPOSCAR"))
        # 备份原始结构
        print("Note: -------------------- ")
        print(f"    finish back up origin inputed structure file into workpath:\n         {self.work_path}")
        if not self.work_path.joinpath(self.input_file_path.name).exists():
            # shutil.copy(self.input_file_path, self.work_path)
            backup_inputfile = self.work_path.joinpath("origin-"+self.input_file_path.name)
            os.system(f"cp -f {self.input_file_path} {backup_inputfile}")

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
        self.cell_parameters = self.get_cell(pstruct) 
        # 获得原子分数坐标
        self.fractional_sites = self.get_coords(pstruct)
        # 获得倒格矢
        # 注意，这是用于固体物理的标准倒数晶格，因数为2π
        self.reciprocal_plattice = self.get_reciprocal_lattice()
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
        self.final_choosed_pp = [] # 最终用来写输入文件用到的赝势都存在这里
        for species in self.species:
            species_name = species.name
            dst_pps = get_pps_for_a_element(species_name, pp_files) # 从已经准备好的pp目录里面拿出来目标赝势
            if len(dst_pps) > 1:
                # 如果目标元素已经有多个赝势放在pp目录中. 例如S元素的三个赝势都放在了pp目录中，你需要手动选择哪个赝势用作计算
                choosed_flag  = False
                while not choosed_flag:
                    for i, p in enumerate(dst_pps):
                        print(f"{i+1}.  ", p)
                    print("\nNote: --------------------")
                    choosed_pp_number = int(input("    please input you want one(Enter the previous number, such as 1,2,3...)\n"))
                    choosed_pp = dst_pps[choosed_pp_number-1]
                    if choosed_pp in dst_pps:
                        choosed_flag = True
                        self.final_choosed_pp.append(choosed_pp)
                    else:
                        choosed_flag = False
                print(self.final_choosed_pp)
            elif len(dst_pps) == 1:
                # 如果目标元素只有1个赝势放在pp目录中. 直接将其加入self.final_choosed_pp即可
                self.final_choosed_pp.append(dst_pps[0])
            elif len(dst_pps) == 0:
                # 如果目标元素在pp目录中一个都没有，就从之前设置的赝势库中选择一个
                print("\nNote: --------------------")
                print(f"    To prepare {species_name} uspp! ")
                dst_pp = self.get_single_uspp(species_name, workpath_pppath)
                self.final_choosed_pp.append(dst_pp)
            else:
                print(f"find many POTCARs {dst_pps}")
                sys.exit(1)

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
            print("    You may not set the potcar_source_libs, you can set it in `qebin.py` ")
            sys.exit(1)

        targetppfiles    = get_pps_for_a_element(species_name, pp_files)
        targetppnames    = [pp for pp in targetppfiles]
        choosed_flag  = False
        while not choosed_flag:
            for i, p in enumerate(targetppnames):
                print(f"{i+1}.  ", p)
            choosed_pp_number = int(input("    please input you want one(Enter the previous number, such as 1,2,3...)\n"))
            choosed_pp = targetppnames[choosed_pp_number-1]
            if choosed_pp in targetppnames:
                src_pp = Path(qe_USPP).joinpath(choosed_pp)
                dst_pp = Path(workpath_pppath).joinpath(choosed_pp)
                shutil.copy(src_pp, dst_pp)
                choosed_flag = True
            else:
                choosed_flag = False
            
        return choosed_pp

    def get_reciprocal_lattice(self):
        """
        读取scffit.out文件获得其中的reciprocal lattice
        Read the scffit.out file to get the reciprocal lattice
        """

        scffit_out_path = self.work_path.joinpath("scffit.out")
        scf_out_path = self.work_path.joinpath("scf.out")
        relax_out_path = self.work_path.joinpath("relax.out")

        if scffit_out_path.exists():
            res_path = scffit_out_path
        elif scf_out_path.exists():
            res_path = scf_out_path
        elif relax_out_path.exists():
            res_path = relax_out_path
        else:
            print("\nNote: --------------------")
            print("    scffit.out, scf.out and relax.out all don't exist. The program can't get reciprocal lattice from them.")
            # time.sleep(2)
            return None

        try:
            # 从scffit.out中获得alat
            alat = float(os.popen(f"sed -n '/lattice parameter (alat)/p' {res_path}").read().split()[4])

            unit_reciprocal_axis = 2*np.pi/alat
            print("\nNote: --------------------")
            print(f"    unit_reciprocal_axis = 2pi/alat=2pi/{alat}={unit_reciprocal_axis}, the unit of alat is `1 a.u.`=0.529117 A")
            print(f"    However !!!!! When you get cartesian coordinates of high symmetry points, unit_reciprocal_axis={1}. Only in this way can we guarantee the consistency of the coordinates of the q points !!!!")
            reciprocal_lattice = []
            for i in range(1,4):
                b = os.popen(f"sed -n '/b({i})/p' {res_path}").read()
                b = re.findall(r"-?\d+\.\d+", b)
                b = [float(bi) for bi in b]  # 这里不需要乘以 unit_reciprocal_axis
                reciprocal_lattice.append(b)
            return np.array(reciprocal_lattice)
        except:
            print("\nNote: --------------------")
            print("    The program fail to get reciprocal lattice from scffit.out, scf.out and relax.out.")
            return None

    def get_cell(self, struct):
        cell_parameters = None
        if self.input_file_path.name == "relax.out" and self.input_file_path.exists():
            awk_order = "awk '/Begin final coordinates/,/End final coordinates/{print $0}'" + f" {self.input_file_path.absolute()} " 
            content = os.popen(awk_order).readlines()
            for idx, line in enumerate(content):
                if "CELL_PARAMETERS (angstrom)" in line:
                    cell_parameters = content[idx+1:idx+4]
                    # 将里面的\n剔除
                    cell_parameters = [cell.strip("\n") for cell in cell_parameters]
            return cell_parameters
        elif self.input_file_path.name == "POSCAR" or self.input_file_path.name == "CONTCAR" or ".vasp" in self.input_file_path.name:
            sed_order = "sed -n '3,5p' " + f" {self.input_file_path.absolute()} "
            content = os.popen(sed_order).readlines()
            cell_parameters = [cell.strip("\n") for cell in content]
            return cell_parameters
        else:
            print("You didn't specify relax.out as inputfile")
            print("So We will get cell-information in the way of PYMATGEN")
            print("cell_parameters")
            cell_parameters = struct.lattice.matrix
            cell_parameters = ['{:>20.16f}    {:>20.16f}    {:>20.16f}'.format(cell[0], cell[1], cell[2]) for cell in cell_parameters]
            return cell_parameters

    def get_coords(self, struct):
        fractional_sites = None
        if self.input_file_path.name == "relax.out" and self.input_file_path.exists():
            awk_order = "awk '/Begin final coordinates/,/End final coordinates/{print $0}'" + f" {os.path.abspath(self.input_file_path)} " 
            content = os.popen(awk_order).readlines()
            for idx, line in enumerate(content):
                if "ATOMIC_POSITIONS (crystal)" in line:
                    fractional_sites = content[idx+1:-1]
                    # 将里面的\n剔除
                    fractional_sites = [coords.strip("\n") for coords in fractional_sites]
            return fractional_sites
        elif self.input_file_path.name == "POSCAR" or self.input_file_path.name == "CONTCAR" or ".vasp" in self.input_file_path.name:
            sed_order1 = "sed -n '6p' " + f" {self.input_file_path.absolute()} "
            elements   = os.popen(sed_order1).read().strip('\n').split()
            sed_order2 = "sed -n '7p' " + f" {self.input_file_path.absolute()} "
            numatoms   = os.popen(sed_order2).read().strip('\n').split()
            # 检查元素的个数与原子种类的个数是否相同
            assert len(elements) == len(numatoms)
            # 检查该POSCAR是否是分数坐标
            sed_order3 = "sed -n '8p' " + f" {self.input_file_path.absolute()} "
            coordstype = os.popen(sed_order3).read().strip('\n')
            assert coordstype[0] == 'D' or coordstype[0] == 'd'
            # 将元素按照原子个数重复写进列表elementlist里
            sed_order4 = "sed -n '9,$p' " + f" {self.input_file_path.absolute()} "
            coordnates = os.popen(sed_order4).readlines()
            element_names= []
            for na, ele in zip(numatoms, elements): # extend() 函数用于在列表末尾一次性追加另一个序列中的多个值
                element_names.extend([ele]*int(na))
            
            # 建立元素列表elementlist与分数坐标的一一对应关系
            fractional_sites = list(zip(element_names, coordnates))
            # 按照self.species的元素顺序去写原子坐标，因为赝势的顺序就是self.species的顺序
            # 这样就可以保证赝势的顺序pp_order和原子坐标的顺序一致！！！！
            pp_order = [spe.name for spe in self.species]
            fractional_sites = sorted(fractional_sites, key=lambda item: pp_order.index(item[0]))
            fractional_sites = [
                '{:<4}   {}'.format(ele, site.strip("\n")) for ele, site in fractional_sites]
            return fractional_sites
        else:
            print("You didn't specify relax.out as inputfile")
            print("So We will get coords-information in the way of PYMATGEN")
            fractional_sites = struct.sites
            element_names = [re.search(r"[A-Za-z]+", str(site.species)).group() for site in fractional_sites]
            fractional_sites = [
                '{:<4}    {:>20.16f}    {:>20.16f}    {:>20.16f}'.format(
                ele,
                site.frac_coords[0], 
                site.frac_coords[1], 
                site.frac_coords[2]) for ele, site in zip(element_names, struct.sites)
                ]
            return fractional_sites
        
def get_pps_for_a_element(species_name:str, pp_files:list):
    """
    species_name: 元素名
    pp_files    : 存储着赝势名的列表
    该函数可以从存储着赝势名的列表中挑选出指定元素的全部赝势
    支持识别的赝势名称格式分别为: 
        H.pbe-van_bm.UPF.txt  首字母大写的元素名.****
        h_pbe_v1.uspp.F.UPF   首字母小写的元素名_****
    """
    # 这里写的实在是太精妙了。我用2个花括号套在一起解决了正则表达式传入花括号的问题
    namelenth = len(species_name)
    if namelenth == 1:
        qepp_name_patter = re.compile(f'''^[{species_name.capitalize()}|{species_name.lower()}][\.|\_]''')
    elif namelenth == 2:

        qepp_name_patter = re.compile(f'''^([{species_name[0].capitalize()}|{species_name[0].lower()}]{species_name[1]})[\.|\_]''')
    else:
        print(f"Inconceivable that the length of the element name is not 1 or 2, it is {namelenth}, its name is {species_name}")
    dst_pps          = list(filter(lambda file: qepp_name_patter.findall(file), pp_files))

    return dst_pps