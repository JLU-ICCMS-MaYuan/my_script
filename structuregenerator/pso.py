import os, sys, shutil
additional_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(additional_path)

import logging
from pathlib import Path
from pprint import pprint
from itertools import combinations
from typing import *

import numpy as np
from ase import Atoms
from ase.geometry.analysis import Analysis
from ase.io import ParseError, read, write
from ase.spacegroup.spacegroup import get_spacegroup
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.core.structure import Structure
from pymatgen.io.ase import AseAtomsAdaptor
from pyxtal.lattice import Lattice
from pyxtal import pyxtal
from pyxtal.tolerance import Tol_matrix

from vasp.vasptools.parser_vasp import VASPOutcarFormat, VASPPoscarFormat
from psolib.utils.convert import dict2Atoms
from psolib.utils.get_enthalpy import get_enthalpy
from psolib.utils.sort_atoms import sort_atoms
from psolib.finger.mixins import FingerPrintMixin, UpdateBestMixin
# from psolib.lib_so.f90sym3dgenerator import Spgcrylat


logger = logging.getLogger(__name__)

def is_bond_reasonable(atoms, type_map, distance_matrix):
    '''
    type_map=['Li','La','H']
    return: min dis of different type bonds including
             `Li-La`, `Li-H`, `La-H`, `Li-Li`, `La-La`, `H-H`
    '''
    number_of_species = distance_matrix.shape[0]
    distance_vector = []
    temp_num = -number_of_species
    for idx in range(number_of_species - 1):
        temp_num += 1
        distance_vector.extend(distance_matrix[idx][temp_num:].tolist())
    distance_vector.extend(np.diagonal(distance_matrix).tolist())
    # distance_vector = [LiLa, LiH, LaH, LiLi, LaLa, HH]

    diff_bond_type = list(map(list, combinations(type_map, 2)))
    for same_type in type_map:
        diff_bond_type.extend([[same_type, same_type]])

    bond_situation = []
    for each_type in diff_bond_type:
        # bond_situation_key = ''.join(each_type)
        ana = Analysis(atoms)
        bond_indice = ana.get_bonds(each_type[0], each_type[1], unique=True)
        if bond_indice == [[]]:
            min_dis = np.array(203)
            # min_dis = np.nan
        else:
            bond_value = ana.get_values(bond_indice)
            min_dis = np.nanmin(bond_value[0])
        bond_situation.append(float(min_dis))

    check_distance = False
    contract_distance = zip(np.array(distance_vector), np.array(bond_situation))
    for standard, current in contract_distance:
        if current <= standard:
            check_distance = False
        else:
            check_distance = True

    return check_distance

class pso(UpdateBestMixin):

    def __init__(
        self,
        work_path,
        popsize,
        distancematrix,
        critic,
        fingerprint,
        nameofatoms,
        threshold,
        lbest,
        maxstep,
        pso_ltype: List[str],
        spacegroup_number: int,
        pso_ratio: float,
        specifywps=None,
        splitwps=None,
        # id_list,
        # pbest_list,
        # gbest,
        # current_atoms__list,
        # fp_mats,
    ):
        self.work_path = work_path

        self.popsize = popsize
        self.distancematrix = np.array(distancematrix)

        self.nameofatoms = nameofatoms
        self.lbest = lbest
        self.maxstep = maxstep
        self.pso_ltype = pso_ltype
        self.spacegroup_number = spacegroup_number
        self.pso_ratio = pso_ratio
        self.specifypwps = specifywps 
        self.splitwps = splitwps
        if self.specifypwps is not None:
            self.group = self.specifypwps._group
        elif self.splitwps is not None:
            self.group = self.splitwps._group
        # self.next_step 表示下一代的'代编号'
        self.next_step = self.get_next_step(work_path)
        # self.current_step 表示当前代的'代编号', 获得它是为了self.generate_step()方法使用
        self.current_step = self.next_step -1 
        self.last_step = self.next_step -1 -1 

        UpdateBestMixin.__init__(self, critic, fingerprint, self.nameofatoms, threshold)

        self.data = [[None, None] for _ in range(popsize)] # 避免浅赋值
        self.fp_mats = (None)
        self.pbest = [[] for _ in range(popsize)]
        self.gbest = {} # gbest 中每个化学配比保存的atoms类数量 由 lbest决定
        self.collect_ini_opt(self.work_path)
        # last_gbest 表示上一代的'全局极小值信息', 例如:
        # gbest_4.extxzy 表示第4代结构经过优化以及综合各方面信息考虑, 
        # 得到的能量最小的结构都存放在这里, 每个化学配比都有一个结构
        last_gbest = Path(self.work_path).joinpath(f"gbest_{self.last_step}.extxyz")
        self.update_current_gbest(last_gbest)
        current_gbest = Path(self.work_path).joinpath(f"gbest_{self.current_step}.extxyz")
        self.store_current_gbest(current_gbest)
        self.store_current_struct(self.work_path)

        self.update_next_step(self.work_path)

        # generate_step 执行时, 使用了self.current_step 变量
        pso_structure = self.generate_step(
            [column for column in range(len(self.data))],
            [self.pbest[column] for column in range(len(self.data))],
            self.gbest,
            [self.data[column] for column in range(len(self.data))],
            self.fp_mats,
        )
        if pso_structure:
            self.store_next_struct(self.work_path, pso_structure)
        else:
            logger.error(f"pso_structure: {pso_structure}")


    @classmethod
    def init_from_config(cls, config_d: dict, *args, **kwargs):

       self = cls(
            work_path=config_d["work_path"],
            popsize=config_d["popsize"],
            distancematrix=config_d["distancematrix"],
            critic=config_d["critic"],
            fingerprint=config_d["fingerprint"],
            nameofatoms=config_d["nameofatoms"],
            threshold=config_d["simthreshold"],
            lbest=config_d["lbest"],
            maxstep=config_d["maxstep"],
            pso_ltype=config_d["pso_ltype"],
            spacegroup_number=config_d["spacegroup_number"],
            pso_ratio=config_d["pso_ratio"],

            specifywps=kwargs.get("specifywps", None),
            splitwps=kwargs.get("splitwps", None),
        #    id_list,
        #    pbest_list, 
        #    gbest, 
        #    current_atoms_list, 
        #    fp_mats
       )
       return self

    def get_next_step(self, work_path):
        '''
        获得下一代结构的编号, 即step值。
        比如:
            1. 产生了第1代结构, 那么对应的next_step应该为2, next_step表示下一代的'代数编号'
               如果没有step这个文件, 那么程序就自己写step文件, 并且将2写进去。
            2. 如果 next_step = 6, 那么说明当前目录下的POSCAR, OUTCAR, CONTCAR都是第5代的结构信息, 下一步我们将向第6代演化
        输入:
            工作目录路径
        返回:
            next_step: 下一代的'代编号'
        '''
        step_file = Path(work_path).joinpath("step")
        if not step_file.exists():
            # raise FileExistsError("step file doesn't exist")
            logger.info("There is no `step` file. So the program will create by itself and write `2` into it!")
            with open(step_file, "w") as f:
                f.write("2")
        
        next_step = eval(open(step_file, "r").read().strip("\n"))

        return next_step

    def collect_ini_opt(self, work_path):
        """
        Function: 
            step_1: In according to the the value of `popsize` in the `pso.ini` file, 
                    the program will determine how many structures will participate in 
                    the PSO structure evolution 
            step_2: In the FOR-cycle, check for the presence of POSCAR and OUTCAR of 
                    every structure.  
                    2.1 the program will convert POSCAR file into atoms_poscar, 
                    2.2 the program will convert OUTCAR file into atoms_outcar, 
                            get the `enthalpy` of every structures, add its `enthalpy` 
                            to the dictionary of `info['enthalpy']` of `Atoms` class.
                            get the `fingerprint` of every structures, 
                            add its `fingerprint` and `fp_type` to the dictionary of `info` of `Atoms` class by blow:
                                self.set_fp(atoms_outcar)
                                the method `set_fp` will add `fp_type` to the dictionary of `info['fp_type']` of `Atoms`
                                the method `set_fp` will add `fingerprint` to the dictionary of `info['fingerprint']` of `Atoms`
                    *** But if OUTCAR is something wrong, the program will use the 
                    corresponding POSCAR to subsitute the OUTCAR. the enthalpy will
                    be set as 610612509 by itself.***
            step_3: set up `self.pbest`. Suppose the POSCAR and OUTCAR of 5 structures are read in. then 
                    the `self.data` will be:
                    self.data = [
                        [atoms_poscar1, atoms_outcar1],
                        [atoms_poscar2, atoms_outcar2],
                        [atoms_poscar3, atoms_outcar3],
                        [atoms_poscar4, atoms_outcar4],
                        [atoms_poscar5, atoms_outcar5],
                    ]
            step_4: set up `self.pbest`. Suppose the POSCAR and OUTCAR of 5 structures are read in. then 
                    the `self.pbest` will be:
                    self.pbest = [
                        [atoms_outcar1],
                        [atoms_outcar2],
                        [atoms_outcar3],
                        [atoms_outcar4],
                        [atoms_outcar5],
                    ]
            step_5: set up `self.fp_mats`
                    current_fp = atoms_outcar.info['fingerprint'], 
                        the `current_fp` is a 3-dim array, its shape is (11, 3, 3)
                    self.fp_mats = np.expand_dims(current_fp, axis=0)
                        the `self.fp_mats` is a 4-dim array, its shape is (x, 11, 3, 3)
                    Suppose the POSCAR and OUTCAR of 5 structures are read in. 
                        then the shape of `self.fp_mats` will be: (5, 11, 3, 3)
                
        input:
            work_path : the parent directory path of poscar_ini, poscar_opt 
        return:
            assign the value for the self.data
        """
        for col in range(self.popsize):
            poscar = Path(work_path).joinpath(f"POSCAR_{col+1}")
            outcar = Path(work_path).joinpath(f"OUTCAR_{col+1}")
            contcar= Path(work_path).joinpath(f"CONTCAR_{col+1}")


            # read the POSCAR and convert it to `atoms_poscar`
            if poscar.exists():
                try:
                    data = VASPPoscarFormat().from_poscar(file_name=poscar)
                    atoms_poscar = dict2Atoms(data, self.nameofatoms)
                except Exception as e:
                    raise ValueError(f"No.{col} structure's POSCAR has problems !!!")
            else:
                raise FileExistsError(f"No.{col+1} poscar didn't exist! So The program can't read it!!!")

            # check whether the OUTCAR file exists or not
            if outcar.exists():
                # if the OUTCAR exists 
                # try to read the OUTCAR and convert it to `atoms_outcar`
                # if read the OUTCAR fails, then the program will read the corresponding POSCAR instead of its OUTCAR
                try:
                    data = VASPOutcarFormat().from_outcar(file_name=outcar)
                    atoms_outcar = dict2Atoms(data, self.nameofatoms)
                    enth = get_enthalpy(outcar, contcar)
                    atoms_outcar.info['enthalpy'] = enth
                    # !!! get_enthalpy函数在读取OUTCAR并返回每原子的焓值时，可能返回nan。
                    # 这是因为有些结构优化起来非常困难，所以优化一半就被杀掉了，所以虽然有该结构对应的OUTCAR
                    # 但它是一个不完成的OUTCAR，读取是可以正常读取的，但是却得不到一个正常的焓值。
                    # 但是我又不希望它最终返回一个nan，因此我设定get_enthalpy读取失败后返回 `610612509`
                    #atoms_outcar.info['sid']      = sid
                    atoms_outcar.info['column']   = int(col)
                    self.set_fp(atoms_outcar)
                except Exception as e:
                    logger.info(f"No.{col} structures's OUTCAR has problems, its POSCAR will be tried to substitue it! ")
                    try:
                        data = VASPPoscarFormat().from_poscar(file_name=poscar)
                        atoms_outcar = dict2Atoms(data, self.nameofatoms)
                        atoms_outcar.info['enthalpy'] = 610612509
                        #atoms_outcar.info['sid']      = sid
                        atoms_outcar.info['column']   = int(col)
                        self.set_fp(atoms_outcar)
                    except Exception as e:
                        raise ValueError(f"No.{col} substituted POSCAR still has problems !!!")
            # if the OUTCAR doesn't exist, then the program will use the corresponding POSCAR instead of its OUTCAR just like beforce.
            else:
                logger.info(f"No.{col+1} structures's OUTCAR doesn't exist, its POSCAR will be tried to substitue it! ")
                try:
                    data = VASPPoscarFormat().from_poscar(file_name=poscar)
                    atoms_outcar = dict2Atoms(data, self.nameofatoms)
                    atoms_outcar.info['enthalpy'] = 610612509
                    #atoms_outcar.info['sid']      = sid
                    atoms_outcar.info['column']   = int(col)
                    self.set_fp(atoms_outcar)
                except Exception as e:
                    raise ValueError(f"No.{col} substituted POSCAR still has problems !!!")

            logger.info(f"finish read No.{col + 1 } structure")

            if isinstance(atoms_poscar, Atoms):
                self.data[col][0] = atoms_poscar
            else:
                raise ValueError(f"atoms_poscar has problem! Perhaps No.{col} POSCAR is something wrong !")
            if isinstance(atoms_outcar, Atoms):
                self.data[col][1] = atoms_outcar
                self.pbest[col].append(atoms_outcar)
            else:
                raise ValueError(f"atoms_outcar has problem! Perhaps No.{col} OUTCAR is something wrong !")

            current_fp = atoms_outcar.info['fingerprint'] # current_fp.shape = (11, 3, 3)
            if self.fp_mats is None:
                self.fp_mats = np.expand_dims(current_fp, axis=0) # self.fp_mats = (1, 11, 3, 3)
            else:
                current_fp = np.expand_dims(current_fp, axis=0)
                self.fp_mats = np.concatenate([self.fp_mats, current_fp])
 
    def update_current_gbest(self, last_step_gbest):
        '''
        更新当前代结构的全局极小值。
        第一步: 读取上一代结构的全局极小值文件: gbest_X.extxyz, 将全部信息存储在self.gbest 字典中
        第二步: 根据当前代结构优化后的结构信息和能量(都存储在self.data中), 更新self.gbest 字典
        输入:
            last_step_gbest: 上一代gbest文件目录
        功能:
            更新self.gbest字典
        '''
        # check the existance of the global minima file of laste step 
        if Path(last_step_gbest).exists():
            logger.info(f"The program will read the {self.last_step}_step gbest and generate the {self.current_step}_step gbest")
            old_gbest_list = read(last_step_gbest, ':', 'extxyz')
            for col, _atoms in enumerate(old_gbest_list):
                _atoms = sort_atoms(atoms=_atoms, elems=self.nameofatoms)
                # try to obtain the `_atoms` named `str(_atoms.symbols)`
                # if there is no the `_atoms` named `str(_atoms.symbols)`, then the program will create an empty list `[]` for this symbols
                # if there is    the `_atoms` named `str(_atoms.symbols)`, then the program will get the List of the corresponding symbols. 
                self.gbest[str(_atoms.symbols)] = self.update_best(
                    self.gbest.get(str(_atoms.symbols), []) + [_atoms], # self.gbest.get(str(_atoms.symbols), []) will return a list which stores the all `_atoms` of this symbols.
                    self.lbest, # `self.lbest` will determine the number of structures about this symbols.
                )
        else:
            logger.info("The program will read the 1_step structure and generate the 2_step structure")

        for col, atoms in enumerate(self.data):
            self.pbest[col] = self.update_best(self.pbest[col] + [atoms[1]], 1)
            self.gbest[str(atoms[1].symbols)] = self.update_best(
                self.gbest.get(str(atoms[1].symbols), []) + [atoms[1]], self.lbest
            )

        logger.info(f"Now there are {len(self.gbest.keys())} stoichiometry")
        logger.info(f"They are respectively:")
        pprint(list(self.gbest.keys()))

    def update_next_step(self, work_path):
        """
        在完成pso结构演化后, self.next_step 指代的那一代初始结构就都产生出来了,
        因此需要更新step文件中的数字, 方便下一次程序产生结构时了解它将产生的是第几代结构
        输入:
            工作目录路径, 该路径下存储着step文件
        功能:
            给self.next_step + 1 => next_next_step, 更新到step文件中, step文件中存储的代数目始终是指向下一代,
            目的就是方便程序读入step后即可知: 下一代是第几代
        """
        step_file = Path(work_path).joinpath("step")
        if not step_file.exists():
            raise FileExistsError("step file doesn't exist")
        with open(step_file, "w") as f:
            next_next_step = self.next_step + 1
            f.write(str(next_next_step))

    def store_current_gbest(self, current_step_gbest):
        for symbols, gbest_list in self.gbest.items():
            write(current_step_gbest, gbest_list, 'extxyz', append=True, parallel=False)
        
    def store_current_struct(self, work_path):
        """
        store the current structures to the directory `result`. For example, 
            `step` file stores number 2, therefore, self.next_step=2, self.current_step=1
            `step` file stores number 5, therefore, self.next_step=5, self.current_step=4
        The execution logic of this instance method:
            1. check the existance of file named `"step_"+str(self.current_step)` , if exist, then remove it.
            2. check the existance of file named `"result` , if not exist, then create it.
            3. extract the enthalpy of the OUTCAR to `atoms_outcar` 
            4. write the file pso_ini_step, pso_opt_step, pso_sor_step
            5. write the file `struct.dat`
                In the `struct.dat`, enthalpy, POSCAR, CONTCAR will be wroten in it by the order mentioned before. 
            6. move all the POSCAR, OUTCAR, CONTCAR into the `current_step_dir`
        """
        # 1. check the existance of file named `"step_"+str(self.current_step)` , if exist, then remove it.
        current_step_dir = Path(work_path).joinpath("step_"+str(self.current_step))
        if current_step_dir.exists():
            os.system(f"rm -fr {str(current_step_dir)}")
        current_step_dir.mkdir()
        # 2. check the existance of file named `"result` , if not exist, then create it.
        result = Path(work_path).joinpath("result")
        if not result.exists():
            os.mkdir(result)

        pso_ini = result.joinpath("pso_ini_" + str(self.current_step))
        pso_opt = result.joinpath("pso_opt_" + str(self.current_step)) 
        pso_sor = result.joinpath("pso_sor_" + str(self.current_step))
        struct_dat = result.joinpath("struct.dat")

        for col in range(self.popsize):
            poscar = Path(work_path).joinpath(f"POSCAR_{col+1}")
            outcar = Path(work_path).joinpath(f"OUTCAR_{col+1}")
            contcar= Path(work_path).joinpath(f"CONTCAR_{col+1}")
            if not poscar.exists() :
                logger.warning(f"POSCAR_{col+1} didn't exist! So the program can't store POSCAR_{col+1}  ")
            if not outcar.exists():
                logger.warning(f"OUTCAR_{col+1} didn't exist! So the program can't store OUTCAR_{col+1}  ")
            if not contcar.exists():
                logger.warning(f"CONTCAR_{col+1} didn't exist! So the program can't store CONTCAR_{col+1} ")

            # 3. extract the enthalpy of the OUTCAR to `atoms_outcar` 
            data = VASPPoscarFormat().from_poscar(file_name=poscar)
            atoms_poscar = dict2Atoms(data, self.nameofatoms)
            try:
                data = VASPOutcarFormat().from_outcar(file_name=outcar)
                atoms_outcar = dict2Atoms(data, self.nameofatoms)
                enth = get_enthalpy(outcar, contcar)
                atoms_outcar.info['enthalpy'] = enth
            except Exception as e:
                data = VASPPoscarFormat().from_poscar(file_name=poscar)
                atoms_outcar = dict2Atoms(data, self.nameofatoms)
                enth = 610612509
                atoms_outcar.info['enthalpy'] = enth

            # 4. write the file pso_ini_step, pso_opt_step, pso_sor_step
            # write pso_ini_step file, such as pso_ini_1,  pso_ini_2,  pso_ini_3, .......
            #   <1> : write `Atoms class` to the file named `"pso_opt_" + str(self.current_step)`
            atoms_poscar.write(pso_ini, format="vasp", append=True)
            # write pso_opt_step file, such as pso_opt_1,  pso_opt_2,  pso_opt_3, .......
            #   <1> write `enthalpy`    to the file named `"pso_opt_" + str(self.current_step)`
            #   <2> : write `Atoms class` to the file named `"pso_opt_" + str(self.current_step)`
            #   For example:
            #       0.00324 
            #       SiO2
            #       1.0
            #       5.000  0.000  0.000
            #       0.000  5.000  0.000
            #       0.000  0.000  5.000
            #       Si O
            #       1  2
            #       0.0    0.0    0.0
            #       0.5    0.5    0.5
            #       0.25   0.25   0.25
            #       610612509 
            #       SiO
            #       1.0
            #       5.000  0.000  0.000
            #       0.000  5.000  0.000
            #       0.000  0.000  5.000
            #       Si O
            #       1  2
            #       0.0    0.0    0.0
            #       0.5    0.5    0.5
            #       .........
            with open(pso_opt, "a") as f:
                f.write(f"{enth}\n")
            atoms_outcar.write(pso_opt, format="vasp", append=True)
            # write pso_opt_step file, such as pso_opt_1,  pso_opt_2,  pso_opt_3, .......
            #   <1> write `enthalpy`    to the file named `"pso_opt_" + str(self.current_step)`
            #   For example:
            #       0.00324 
            #       610612509 
            with open(pso_sor, "a") as f:
                f.write(f"{enth}\n")

            # 5. write the file `struct.dat`
            #   In the `struct.dat`, enthalpy, POSCAR, CONTCAR will be wroten by the order mentioned before.
            #   write poscar information
            with open(struct_dat, "a") as f:
                nstruct = (self.current_step-1)*self.popsize + (col + 1) # add another `1` due to python beginning from 0
                f.write("        nstruct= {:<10}\n".format(nstruct))
                f.write("++++++++++++++++++++++++++++++++++\n")
                f.write("        Initial Structure\n")
                spggroup_symbol = get_spacegroup(atoms_poscar).symbol
                f.write("Space Group= {:<10}\n".format(spggroup_symbol))
                volume = atoms_poscar.get_volume()
                f.write("Volume= {:<.8f}\n".format(volume))
                f.write("Number Species= {:<20}\n".format(str(len(self.nameofatoms))))
                symbols = atoms_poscar.symbols
                _numions =[symbols.count(ele) for ele in self.nameofatoms]
                numions = list(map(str, _numions))
                f.write("Ele_Num= {:<50}\n".format(' '.join(numions)))
                f.write(" ---------------------------------\n")
                f.write("lat_matrix\n")
                for vector in atoms_poscar.cell:
                    f.write("{:<.8f} {:<.8f} {:<.8f}\n".format(vector[0], vector[1], vector[2]))
                f.write(" ---------------------------------\n")
                f.write(" ---------------------------------\n")
                f.write("Atomic Positions\n")
                fractional_position = atoms_poscar.get_scaled_positions()
                for vector in fractional_position:
                    f.write("{:<.8f} {:<.8f} {:<.8f}\n".format(vector[0], vector[1], vector[2]))
                f.write(" ---------------------------------\n")
            #   write outcar information
            with open(struct_dat, "a") as f:
                f.write(" ======================================\n")
                f.write("        Optimized Structure\n")
                f.write("Energy= {:<.10f}\n".format(atoms_outcar.info['enthalpy']))
                volume = atoms_outcar.get_volume()
                f.write("Volume= {:<.8f}\n".format(volume))
                f.write("Number Species= {:<20}\n".format(str(len(self.nameofatoms))))
                symbols = atoms_outcar.symbols
                _numions =[symbols.count(ele) for ele in self.nameofatoms]
                numions = list(map(str, _numions))
                f.write("Ele_Num= {:<50}\n".format(' '.join(numions)))
                f.write(" ---------------------------------\n")
                f.write("lat_matrix\n")
                for vector in atoms_outcar.cell:
                    f.write("{:<.8f} {:<.8f} {:<.8f}\n".format(vector[0], vector[1], vector[2]))
                f.write(" ---------------------------------\n")
                f.write(" ---------------------------------\n")
                f.write("Atomic Positions\n")
                fractional_position = atoms_outcar.get_scaled_positions()
                for vector in fractional_position:
                    f.write("{:<.8f} {:<.8f} {:<.8f}\n".format(vector[0], vector[1], vector[2]))
                f.write(" ---------------------------------\n")
                f.write("\n")
                f.write("\n")
            # 6. move all the POSCAR, OUTCAR, CONTCAR into the `current_step_dir`
            # For No.col structure, its POSCAR, OUTCAR, CONTCAR has been extracted completely !!!
            # Now the program will move its POSCAR OUTCAR CONTCAR into the directory `current_step_dir` !!!
            if poscar.exists(): shutil.move(poscar,  current_step_dir); 
            if outcar.exists(): shutil.move(outcar,  current_step_dir); 
            if contcar.exists(): shutil.move(contcar, current_step_dir); 

    def store_next_struct(self, work_path, atoms_list):
        """
        store the next structures to the directory `work_path`. For example, 
            POSCAR_1, POSCAR_2 .... POSCAR_col .... 
        The execution logic of this instance method:
            write all structures created by PSO method to the `work_path` directory
        """
        for col, atoms in enumerate(atoms_list):
            poscar = Path(work_path).joinpath(f"POSCAR_{col+1}")
            if atoms:
                write(poscar, atoms, format="vasp")
                logger.info(f"try successfully write POSCAR_{col+1}")
            else:
                logger.warning(f"The `Class Atoms` for POSCAR_{col+1} doesn't exist !!! So the program can't write it !!!")


    def generate_step(
        self, 
        id_list: list[int],      # sid   401,    402,    403  ...  全局指标：每代100个结构，401为第4代第1个结构, 402为第4代第2个结构, ...
        pbest_list: list[Atoms], # pbest atoms1, atoms2, atoms...  局域指标：某一代结构的第一个结构, 第二个结构... 
        gbest: dict[Atoms], 
        current_atoms_list: list[Atoms], 
        fp_mats,
        ):

        gen_structures_list = []
        for column, current_atoms in enumerate(current_atoms_list):
            if pbest_list[column]:
                logger.info(f"generate No.{column+1} {pbest_list[column][0].symbols}")
                gen_one = self.generate_one(
                    id_list[column],
                    pbest_list[column],
                    gbest[str(current_atoms[0].symbols)],
                    current_atoms,
                    fp_mats,
                )
                gen_structures_list.append(gen_one)
            else:
                logger.warning(f"pbest_list {pbest_list[column]}")
        
        return gen_structures_list

    def generate_one(
        self, 
        sid, 
        pbest, 
        gbest: list, 
        current_atoms: Atoms, 
        fp_mats
        ) -> Atoms:
        
        random_ratio = np.random.uniform(0,1)
        if random_ratio < self.pso_ratio:
            gen_atoms = self.__pso_gen(pbest, gbest, current_atoms, fp_mats)
        else:
            # gen_atoms = self.__random_gen_method2(current_atoms)
            gen_atoms = self.__random_gen()
        return gen_atoms

    def __pso_gen(
        self,
        pbest: list[Atoms],
        gbest: list[Atoms],
        current_atoms: list[Atoms],
        fp_mats,
    ):
        logger.info(f"create the structure by PSO for {current_atoms[0].symbols}")
        distance_of_ion = self.distancematrix
        init_atoms = current_atoms[0]
        opt_atoms = current_atoms[1]

        omega_max, omega_min = 0.9, 0.4
        iter_max  = self.maxstep
        step = self.current_step + 1  # 更新step
        omega = omega_max - (omega_max - omega_min) * step / iter_max
        c1 = c2 = 2

        dists, indexes, _ = self.get_similarity(
            opt_atoms.info['fingerprint'],
            np.array([atoms.info['fingerprint'] for atoms in gbest]),
        )
        gbest = gbest[indexes[1]]
        pbest_pos = pbest[0].get_scaled_positions()
        # todo: decide which gbest should be used
        gbest_pos        = gbest.get_scaled_positions()
        init_pos         = init_atoms.get_scaled_positions()
        opt_volume       = opt_atoms.get_volume(); opt_symbol = opt_atoms.symbols
        detla_pbest_init = pbest_pos - init_pos
        detla_gbest_init = gbest_pos - init_pos
        for _ in range(5):
            r1 = np.random.rand()
            r2 = np.random.rand()
            init_velocity = init_atoms.arrays.get('caly_velocity', np.random.uniform(0, 1, (len(init_atoms), 3)))
            # np.random.uniform(0, 1, (len(init_atoms), 3)) represents creating a 2-dim array, its shape is len(init_atoms)*3
            velocity = (
                omega * init_velocity
                + c1 * r1 * detla_pbest_init
                + c2 * r2 * detla_gbest_init
            )
            gen_pso_pos = init_pos + velocity # gen_pso_pos is the new coordinates created by PSO method !!!

            # If the range of crystal systems is given in `pso.ini`,
            # the program will randomly select one as `_ltype`, 
            # and then the corresponding lattice matrix will be created. 
            _ltype = np.random.choice(self.pso_ltype)
            _volume = opt_atoms.cell.volume + 20
            lattice = Lattice(ltype=_ltype, volume=_volume).generate_matrix()

            # lattice = np.ones((3, 3), dtype=np.float64, order='F')
            # spacegroupid = np.random.uniform(1, 231, 1)
            # Spgcrylat().spggenlat(opt_volume, spacegroupid, lattice)
            #   opt_volume is the volume of optimized structure
            #   spacegroupid is the space group ID
            #   after Spgcryla().spggenlat() is running, the `lattice` will be assigned new volume !!!
            # logger.warning("This lattice may not be physical.")

            new_atoms = Atoms(
                symbols=opt_symbol,
                scaled_positions=gen_pso_pos,
                cell=lattice,
                pbc=True,
            )
            new_atoms.arrays['caly_velocity'] = velocity
            fp = self.get_fp(new_atoms)
            if (
                is_bond_reasonable(new_atoms, self.nameofatoms, distance_of_ion)
                and not self.get_similarity(fp, fp_mats)[2]  # not is_sim
            ):
                logger.info(f"No.{_+1} attempt for {opt_symbol} succeeded !!!")
                break
            else:
                logger.info(f"No.{_+1} attempt for {opt_symbol} failed !!!")
        else:
            # 
            __ini_atoms = sort_atoms(init_atoms, self.nameofatoms)
            logger.warning(f"The structure may not be physical.")
            logger.warning("But the program will still adopt this structure !!!")

        self.set_fp(new_atoms)
        return new_atoms

    def __same_symbol_gen(self, current_atoms):
        """
        Function:
            create a new structure whose symbols is the same as `current_atoms`
            And do not require the `occupied wyckoff positions` of new_structure are the same as current_atoms[0]
            Additional Notes: 
                current_atoms[0].symbols == atoms_poscar.symbols
        input parameter:
            current_atoms: List[Atoms]. It's a list which dimension is 1 and which size is 2. 
            current_atoms = [atoms_poscar, atoms_outcar]
        """
        # current_atoms is a list whose size is 2
        # current_atoms[0] is the initial   structure whose type is `Atoms`
        # current_atoms[1] is the optimized structure whose type is `Atoms`
        symbols = current_atoms[0].symbols
        species_radius = list(zip(self.nameofatoms, self.distancematrix.diagonal()/2.0))
        tm = Tol_matrix()
        for ele_r1, ele_r2 in combinations(species_radius, 2):
            tm.set_tol(ele_r1[0], ele_r2[0], ele_r1[1]+ele_r2[1])

        struct = pyxtal()
        logger.info(f"create the structure by RANDOM for {current_atoms[0].symbols}")
        for i in range(1000):
            try:
                spacegroup_number = self.random_get_spacegroup(self.spacegroup_number)
                struct.from_random(
                    3,
                    spacegroup_number,
                    species=self.nameofatoms,
                    numIons=[symbols.count(ele) for ele in self.nameofatoms],
                    factor=1.5,
                    tm=tm,
                )
                _new_atoms = struct.to_ase()
                new_atoms = sort_atoms(_new_atoms, self.nameofatoms)
                return new_atoms
            except Exception as e:
                logger.warning(f"create structure No.{i+1} failed !!!")
        else:
            logger.warning(f"Its input information can't create the random {current_atoms[0].symbols}")
            return None

    def __random_gen(self):
        '''
        Function:
            This method of generating structures adopts `specify_wyckoffs._gen_specify_symbols` in `specify_wyckoffs.py` file
            The program will get the `symbols` from `current_atoms[0]`. Then this instance method will create a structure whose 
            stoichiometry and occupied wyckoff positions are the same as current_atoms[0].
            Additional Notes: 
                current_atoms[0].symbols == atoms_poscar.symbols
        input parameter:
            current_atoms: List[Atoms]. It's a list which dimension is 1 and which size is 2. 
            current_atoms = [atoms_poscar, atoms_outcar]
        '''
        if self.specifypwps:
            
            while True:
                stru = self.specifypwps._gen_randomly()
                if hasattr(stru, "is_ordered"): # 判断结构是否是分数占据的无序结构
                    pstru = SpacegroupAnalyzer(stru).get_conventional_standard_structure()
                    _new_atoms = AseAtomsAdaptor.get_atoms(pstru)
                    new_atoms = sort_atoms(_new_atoms, self.nameofatoms)
                    logger.info(f"The program successfully created a structuere by `__random_gen` {new_atoms.symbols}")
                    return new_atoms
        elif self.splitwps:
            while True:
                stru = self.splitwps._gen_randomly()
                if hasattr(stru, "is_ordered"): # 判断结构是否是分数占据的无序结构
                    pstru = SpacegroupAnalyzer(stru).get_conventional_standard_structure()
                    _new_atoms = AseAtomsAdaptor.get_atoms(pstru)
                    new_atoms = sort_atoms(_new_atoms, self.nameofatoms)
                    logger.info(f"The program successfully created a structuere by `__random_gen` {new_atoms.symbols}")
                    return new_atoms

    def random_get_spacegroup(self, spacegroup_number):
        """
        根据给定的空间群所在的晶系, 随机给出一个该晶系下的任意一个空间群
        """
        triclinic    = [i for i in range(1, 3)]
        monoclinic   = [i for i in range(3, 16)]
        orthorhombic = [i for i in range(16, 75)]
        tetragonal   = [i for i in range(75, 143)]
        trigonal     = [i for i in range(143, 168)]
        hexagonal    = [i for i in range(168, 195)]
        cubic        = [i for i in range(195, 231)]

        if spacegroup_number in triclinic:
            select_crystal = triclinic.remove(spacegroup_number)
        elif spacegroup_number in monoclinic:
            select_crystal = monoclinic.remove(spacegroup_number)
        elif spacegroup_number in orthorhombic:
            select_crystal = orthorhombic.remove(spacegroup_number)
        elif spacegroup_number in tetragonal:
            select_crystal = tetragonal.remove(spacegroup_number)
        elif spacegroup_number in trigonal:
            select_crystal = trigonal.remove(spacegroup_number)
        elif spacegroup_number in hexagonal:
            select_crystal = hexagonal.remove(spacegroup_number)
        elif spacegroup_number in cubic:
            select_crystal = cubic.copy()
            select_crystal.remove(spacegroup_number)
        else:
            select_crystal = [i for i in range(1, 231)]

        new_spacegroup_number = np.random.choice(select_crystal)
        
        return new_spacegroup_number