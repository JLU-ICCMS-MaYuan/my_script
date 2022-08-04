'''
qeSuperconductTc.py -pos scripts_tests/POSCAR -caldir scripts_tests/out
-pos    scripts_tests/POSCAR 
-caldir scripts_tests/out
'''

import os
import re
import shutil
import logging
from pathlib import Path
from itertools import chain

from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.core.periodic_table import Element
from pymatgen.io.vasp import Poscar
from pymatgen.io.ase import AseAtomsAdaptor
from ase.io import read


logging.basicConfig(
    level = logging.INFO, 
    format='%(asctime)s | %(name)s | %(levelname)s ---- %(message)s'
    )
logger = logging.getLogger(__name__)

class qe_inputpara:

    def __init__(
        self,
        input_file_path: str,
        work_underpressure: str,
        workpath_pppath: str,
        pressure: int,
        kpoints_dense: list[int],
        kpoints_sparse: list[int],
        qpoints: list[int],
        inserted_points_num: int,
        run_mode: str,
        ):
        self.input_file_path      = input_file_path
        self.work_underpressure   = work_underpressure
        self.workpath_pppath      = workpath_pppath
        self.pressure             = pressure
        self.k1_dense , self.k2_dense , self.k3_dense  = kpoints_dense
        self.k1_sparse, self.k2_sparse, self.k3_sparse = kpoints_sparse
        self.q1,        self.q2,        self.q3        = qpoints
        self.q_total_amount           = None
        self.q_non_irreducible_amount = None
        self.inserted_points_num      = None 
        self.q_coordinate_list        = []
        self.q_weight_list            = [] 
        self.run_mode                 = run_mode
        self.relax_ase  = read(self.input_file_path)
        self.struct     = AseAtomsAdaptor.get_structure(self.relax_ase)
        self.get_struct_info(self.struct)
        # prepare the uspp file in pp directory
        self.get_USPP(self.workpath_pppath)
        # prepare q_total_amount, q_non_irreducible_amount, q_coordinate_list 

        dyn0_names = list(Path(self.work_underpressure).glob("*.dyn0"))
        if len(dyn0_names)==1:
            dyn0_path = str(dyn0_names[0].absolute())
        else:
            logger.warning("No exist *.dyn0! ")

        if self.run_mode == "merge":
            self.get_q_from_dyn0(dyn0_path)
            self.merge(self.work_underpressure)
        if self.run_mode == "lambda":
            q2r_out = list(Path(self.work_underpressure).glob("q2r.out"))
            if q2r_out:
                q2r_path = str(q2r_out[0].absolute()) 
                self.get_q_from_dyn0(dyn0_path, q2r_path=q2r_path)
        if self.run_mode == "matdyn":
            self.inserted_points_num = inserted_points_num
            self.path_name_coords    = self.get_hspp()

    def get_struct_info(self, struct):
        
        spa = SpacegroupAnalyzer(struct)
        # bstruct = spa.get_conventional_standard_structure()
        pstruct = spa.get_primitive_standard_structure()
        Poscar(pstruct).write_file("PPOSCAR")

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
        pp_files = os.listdir(workpath_pppath)
        if len(pp_files) != len(self.species):
            qe_USPP1 = os.path.abspath("/work/home/mayuan/POT/qe-pp/all_pbe_UPF_v1.5")
            qe_USPP2 = os.path.abspath("/public/home/mayuan/POT/qe-pp/all_pbe_UPF_v1.5")
            if os.path.exists(qe_USPP1):
                qe_USPP = qe_USPP1
                ppfiles = os.listdir(qe_USPP1)
            elif os.path.exists(qe_USPP2):
                qe_USPP = qe_USPP2
                ppfiles = os.listdir(qe_USPP2)
            self.final_choosed_pp = []
            for species in self.species:
                species_name = species.name.lower()
                logger.info(f"species_name {species_name}")
                targetppfiles = filter(lambda file: re.search(species_name+"\_", file.lower()), ppfiles)
                targetppnames = [pp for pp in targetppfiles]
                choosed_flag  = False
                while not choosed_flag:
                    choosed_pp = input(f"{targetppnames}, \nplease input you want one\n")
                    if choosed_pp in targetppnames:
                        src_pp = os.path.join(qe_USPP,        choosed_pp)
                        dst_pp = os.path.join(workpath_pppath, choosed_pp)
                        shutil.copy(src_pp, dst_pp)
                        choosed_flag = True
                        self.final_choosed_pp.append(choosed_pp)
                    else:
                        choosed_flag = False
        else:
            self.final_choosed_pp = pp_files
        logger.info(f"the choosed pp is {self.final_choosed_pp}")

    def checkfile(self):
        if "pp" not in os.listdir(self.work_path):
            logger.warning(f"please prepare pp directory!")
        filesordirs = os.listdir(self.work_underpressure)
        if filesordirs:
            if "scf.fit.in" not in filesordirs:
                logger.warning(f"please prepare the scf.fit.in")
            if "scf.in" not in filesordirs:
                logger.warning(f"please prepare the scf.in")
            if "no_split_ph.in" not in filesordirs:
                logger.warning(f"please prepare the no_split_ph.in")
            if "q2r.in" not in filesordirs:
                logger.warning(f"please prepare the q2r.in")
            if "matdyn.in" not in filesordirs:
                logger.warning(f"please prepare the matdyn.in")
            if "matdyn.dos.in" not in filesordirs:
                logger.warning(f"please prepare the matdyn.dos.in")
            if "lambda.in" not in filesordirs:
                logger.warning(f"please prepare the lambda.in")
        else:
            logger.warning("There is no any .in file!!!")

    def get_q_from_scfout(self, dir):
        if not os.path.exists(dir):
            raise FileExistsError ("scf.out didn't exist!")
        content = open(dir, "r").readlines()
        def find_k(item):
            if re.search(r"k\(\s*\d+\)\s*=\s*", item):
                return item

        result = filter(find_k, content)
        for res in result:
            ks  = re.findall(r"\-?\d+\.\d+", res.split(",")[0])
            self.q_coordinate_list.append(ks)
            wp = re.findall(r"\-?\d+\.\d+", res.split(",")[1])
            nqs = float(wp[0]) * self.q_total_amount / 2
            self.q_weight_list.append(nqs)

        self.q_total_amount           = self.q1 * self.q2 * self.q3
        self.q_non_irreducible_amount = len(self.q_coordinate_list)

        return self.q_total_amount,    self.q_non_irreducible_amount, \
               self.q_coordinate_list, self.q_weight_list

    def get_q_from_dyn0(self, dir, q2r_path=None):
        if not os.path.exists(dir):
            raise FileExistsError ("dyn0 file doesn't exist!")
        content = open(dir, "r").readlines()
        _q_total_amount = content[0].strip("\n").split()
        q_total_amount  = list(map(int, _q_total_amount))
        if q_total_amount != [self.q1, self.q2, self.q3]:
            raise ValueError ("q points set wrong")

        def find_q(item):
            if re.search(r"E[\+|\-]", item):
                return item
        def find_q_weight(item):
            if re.search(r"nqs\=", item):
                return item

        self.q_total_amount           = self.q1 * self.q2 * self.q3
        self.q_non_irreducible_amount = content[1]
        _q_coordinate_list            = list(filter(find_q, content))
        self.q_coordinate_list        = [q_string.strip("\n").split() for q_string in _q_coordinate_list]
        # get q_weigth_list
        if q2r_path is not None:
            q2r_out            = open(q2r_path, "r").readlines()
            _q_weight_list     = list(filter(find_q_weight, q2r_out))
            self.q_weight_list = [re.search(r"\d+", qwt).group() for qwt in _q_weight_list]
        
        return  self.q_total_amount,    self.q_non_irreducible_amount, \
                self.q_coordinate_list, self.q_weight_list

    def merge(self, dir):
        elph_dir_path = os.path.join(dir, "elph_dir")
        if not os.path.exists(elph_dir_path):
            os.makedirs(elph_dir_path)
        
        for i in range(int(self.q_non_irreducible_amount)):
            src_elph   = os.path.join(dir, str(i+1), "elph_dir", "elph.inp_lambda.1")
            dst_elph   = os.path.join(elph_dir_path, "elph.inp_lambda."+str(i+1))
            shutil.copy(src_elph, dst_elph)
            logger.info(f"elph.inp_lambda.1 copy finished \n {dst_elph}")

            src_dyn    = os.path.join(dir, str(i+1), self.system_name+".dyn")
            dst_dyn    = os.path.join(dir,           self.system_name+".dyn"+str(i+1))
            shutil.copy(src_dyn, dst_dyn)
            logger.info(f"{self.system_name}.dyn copy finished \n {dst_dyn}")

            for j in range(51, 61):
                src_a2Fq2r = os.path.join(dir, str(i+1), "elph_dir", "a2Fq2r."+str(j)+".1")
                dst_a2Fq2r = os.path.join(elph_dir_path,             "a2Fq2r."+str(j)+"."+str(i+1))
                shutil.copy(src_a2Fq2r, dst_a2Fq2r)
                logger.info(f"a2Fq2r.{str(j)}.1 copy finished \n {dst_dyn}")

    def get_hspp(self):
        """
        This method is to get high symmetry paths and points
        """ 
        lat             = self.relax_ase.cell.get_bravais_lattice()
        path_name_string      = lat.special_path

        _path_name_list = [[ p for p in pp] for pp in path_name_string.split(",")]
        logger.info(f"the high symmetry points path is {_path_name_list}")
        high_symmetry_type = input("please input the mode you want: 'all_points', 'main_points', or Nothing to input\n\n")
        if not high_symmetry_type:
            high_symmetry_type = "all_points" # default
        if "," in path_name_string:
            if high_symmetry_type == "all_points":
                path_name_list = list(chain.from_iterable(_path_name_list))
                logger.info(f"the choosed high symmetry points path is \n {path_name_list}")
            elif high_symmetry_type == "main_points":
                path_name_list = max(_path_name_list, key=len)
                logger.info(f"the choosed high symmetry points path is \n {path_name_list}")
        else:
            path_name_list = [ pp for pp in path_name_string]
            logger.info(f"the choosed high symmetry points path is \n {path_name_list}")

        special_points       = lat.get_special_points()
        path_coords          = [list(special_points[point_name]) for point_name in path_name_list]
        path_name_coords= list(zip(path_name_list, path_coords))
        return path_name_coords 