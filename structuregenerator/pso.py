from json.encoder import py_encode_basestring
from multiprocessing.spawn import old_main_modules
import os, sys, shutil
additional_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(additional_path)

import logging
from pathlib import Path
from pprint import pprint
from itertools import combinations

import numpy as np
from ase import Atoms
from ase.geometry.analysis import Analysis
from ase.io import ParseError, read, write

from vasp.vasptools.parser_vasp import VASPOutcarFormat, VASPPoscarFormat
from psolib.utils.convert import dict2Atoms
from psolib.utils.get_enthalpy import get_enthalpy
from psolib.finger.mixins import FingerPrintMixin, UpdateBestMixin
from psolib.lib_so.f90sym3dgenerator import Spgcrylat

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
        # old_gbest 表示上一代的'全局极小值信息', 例如:
        # gbest_4.extxzy 表示第4代结构经过优化以及综合各方面信息考虑, 
        # 得到的能量最小的结构都存放在这里, 每个化学配比都有一个结构
        last_gbest = Path(self.work_path).joinpath(f"gbest_{self.last_step}.extxyz")
        self.update_current_gbest(last_gbest)
        # generate_step 执行时, 使用了self.current_step 变量
        pso_structure = self.generate_step(
            [column for column in range(len(self.data))],
            [self.pbest[column] for column in range(len(self.data))],
            self.gbest,
            [self.data[column] for column in range(len(self.data))],
            self.fp_mats,
        )

        current_gbest = Path(self.work_path).joinpath(f"gbest_{self.current_step}.extxyz")
        self.update_next_step(self.work_path)
        self.store_current_gbest(current_gbest)
        self.store_current_struct(self.work_path)
        self.store_next_struct(self.work_path, pso_structure)


    @classmethod
    def init_from_config(cls, config_d: dict):
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
            2. 如果 next_step = 6, 那么说明当前目录下的POSCAR, OUTCAR, CONTCAR都是第5代的结构信息, 下一步我们将向第6代演化
        输入:
            工作目录路径
        返回:
            next_step: 下一代的'代编号'
        '''
        step_file = Path(work_path).joinpath("step")
        if not step_file.exists():
            raise FileExistsError("step file doesn't exist")

        next_step = eval(open(step_file, "r").read().strip("\n"))
        return next_step

    def collect_ini_opt(self, work_path):
        """
        input:
            work_path : the parent directory path of poscar_ini, poscar_opt 
        return:
            assign the value for the self.data
        """
        current_symbols = []
        for col in range(self.popsize):
            poscar = Path(work_path).joinpath(f"POSCAR_{col+1}")
            outcar = Path(work_path).joinpath(f"OUTCAR_{col+1}")
            contcar= Path(work_path).joinpath(f"CONTCAR_{col+1}")

            if not poscar.exists() or not outcar.exists():
                logger.warning(f"No.{col} structure didn't exist! ")

            # 给poscar的 atoms类赋值
            data = VASPPoscarFormat().from_poscar(file_name=poscar)
            atoms_poscar = dict2Atoms(data, self.nameofatoms)

            # 给outcar的 atoms类赋值
            try:
                data = VASPOutcarFormat().from_outcar(file_name=outcar)
                atoms_outcar = dict2Atoms(data, self.nameofatoms)
                enth = get_enthalpy(outcar, contcar)
                atoms_outcar.info['enthalpy'] = enth
                #atoms_outcar.info['sid']      = sid
                atoms_outcar.info['column']   = int(col)
                self.set_fp(atoms_outcar)
            except Exception as e:
                data = VASPPoscarFormat().from_poscar(file_name=poscar)
                atoms_outcar = dict2Atoms(data, self.nameofatoms)
                atoms_outcar.info['enthalpy'] = 610612509
                #atoms_outcar.info['sid']      = sid
                atoms_outcar.info['column']   = int(col)
                self.set_fp(atoms_outcar)
            logger.info(f"finish read No.{col + 1 } structure")

            
            self.data[col][0] = atoms_poscar
            self.data[col][1] = atoms_outcar
            self.pbest[col].append(atoms_outcar)

            current_fp = atoms_outcar.info['fingerprint']
            if self.fp_mats is None:
                self.fp_mats = np.expand_dims(current_fp, axis=0)
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
        if Path(last_step_gbest).exists():
            logger.info(f"The program will read the {self.last_step}_step gbest and generate the {self.current_step}_step gbest")
            old_gbest_list = read(last_step_gbest, ':', 'extxyz')
            for col, atoms in enumerate(old_gbest_list):
                self.gbest[atoms.get_chemical_formula('metal')] = self.update_best(
                    self.gbest.get(atoms.get_chemical_formula('metal'), []) + [atoms], self.lbest
                )
        else:
            logger.info("The program will read the 1_step structure and generate the 2_step structure")

        for col, atoms in enumerate(self.data):
            self.pbest[col] = self.update_best(self.pbest[col] + [atoms[1]], 1)
            self.gbest[atoms[1].get_chemical_formula('metal')] = self.update_best(
                self.gbest.get(atoms[1].get_chemical_formula('metal'), []) + [atoms[1]], self.lbest
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
        current_step_dir = Path(work_path).joinpath("step_"+str(self.current_step))
        result = Path(work_path).joinpath("result")
        if current_step_dir.exists():
            os.system(f"rm -fr {str(current_step_dir)}")
        current_step_dir.mkdir()
        if not result.exists():
            os.mkdir(result)

        pso_ini = result.joinpath("pso_ini_" + str(self.current_step))
        pso_opt = result.joinpath("pso_opt_" + str(self.current_step)) 
        pso_sor = result.joinpath("pso_sor_" + str(self.current_step))
        # 获得焓值
        current_symbols = []
        for col in range(self.popsize):
            poscar = Path(work_path).joinpath(f"POSCAR_{col+1}")
            outcar = Path(work_path).joinpath(f"OUTCAR_{col+1}")
            contcar= Path(work_path).joinpath(f"CONTCAR_{col+1}")

            if not poscar.exists() or not outcar.exists():
                logger.warning(f"No.{col} structure didn't exist! ")

            # 给poscar的 atoms类赋值
            data = VASPPoscarFormat().from_poscar(file_name=poscar)
            atoms_poscar = dict2Atoms(data, self.nameofatoms)

            try:
                data = VASPOutcarFormat().from_outcar(file_name=outcar)
                atoms_outcar = dict2Atoms(data, self.nameofatoms)
                enth = get_enthalpy(outcar, contcar)

            except Exception as e:
                data = VASPPoscarFormat().from_poscar(file_name=poscar)
                atoms_outcar = dict2Atoms(data, self.nameofatoms)
                enth = 610612509

            atoms_poscar.write(pso_ini, format="vasp", append=True)
            with open(pso_opt, "a") as f:
                f.write(f"{enth}\n")
            atoms_outcar.write(pso_opt, format="vasp", append=True)
            with open(pso_sor, "a") as f:
                f.write(f"{enth}\n")

            shutil.move(poscar,  current_step_dir); 
            shutil.move(outcar,  current_step_dir); 
            shutil.move(contcar, current_step_dir); 

    def store_next_struct(self, work_path, atoms_list):
        for col, atoms in enumerate(atoms_list):
            poscar = Path(work_path).joinpath(f"POSCAR_{col+1}")
            write(poscar, atoms, format="vasp")


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
            gen_one = self.generate_one(
                id_list[column],
                pbest_list[column],
                gbest[current_atoms[0].get_chemical_formula('metal')],
                current_atoms,
                fp_mats,
            )
            gen_structures_list.append(gen_one)
        return gen_structures_list

    def generate_one(
        self, 
        sid, 
        pbest, 
        gbest: list, 
        current_atoms: Atoms, 
        fp_mats
        ) -> Atoms:
        
        gen_atoms = self.__pso_gen(pbest, gbest, current_atoms, fp_mats)

        return gen_atoms

    def __pso_gen(
        self,
        pbest: list[Atoms],
        gbest: list[Atoms],
        current_atoms: list[Atoms],
        fp_mats,
    ):

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
        opt_volume       = opt_atoms.get_volume()
        detla_pbest_init = pbest_pos - init_pos
        detla_gbest_init = gbest_pos - init_pos
        for _ in range(5):
            print(_)
            r1 = np.random.rand()
            r2 = np.random.rand()
            init_velocity = init_atoms.arrays.get('caly_velocity', np.random.uniform(0, 1, (len(init_atoms), 3)))
            velocity = (
                omega * init_velocity
                + c1 * r1 * detla_pbest_init
                + c2 * r2 * detla_gbest_init
            )
            gen_pso_pos = init_pos + velocity

            spacegroupid = np.random.uniform(1, 231, 1)
            lattice = np.ones((3, 3), dtype=np.float64, order='F')
            Spgcrylat().spggenlat(opt_volume, spacegroupid, lattice)
            # logger.warning("This lattice may not be physical.")

            new_atoms = Atoms(
                symbols=opt_atoms.symbols,
                scaled_positions=gen_pso_pos,
                cell=lattice,
                pbc=True,
            )
            new_atoms.arrays['caly_velocity'] = velocity
            fp = self.get_fp(new_atoms)
            name_of_atoms = self.nameofatoms
            if (
                is_bond_reasonable(new_atoms, name_of_atoms, distance_of_ion)
                and not self.get_similarity(fp, fp_mats)[2]  # not is_sim
            ):
                break
        else:
            logger.warning("This structure may not be physical.")

        self.set_fp(new_atoms)
        return new_atoms