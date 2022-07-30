#!/public/home/mayuan/miniconda3/envs/cage/bin/python3
#!/work/home/mayuan/miniconda3/envs/cage/bin/python3

'''
qeSuperconductTc.py -pos scripts_tests/POSCAR -caldir scripts_tests/out
-pos    scripts_tests/POSCAR 
-caldir scripts_tests/out
'''

import os
import re
import shutil
import logging

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
        ):
        self.input_file_path      = input_file_path
        self.work_underpressure   = work_underpressure
        self.workpath_pppath      = workpath_pppath
        self.pressure             = pressure
        self.k1_dense , self.k2_dense , self.k3_dense  = kpoints_dense
        self.k1_sparse, self.k2_sparse, self.k3_sparse = kpoints_sparse
        self.q1,        self.q2,        self.q3        = qpoints
        relax_ase   = read(self.input_file_path)
        self.struct = AseAtomsAdaptor.get_structure(relax_ase)
        self.get_struct_info(self.struct)
        # prepare the uspp file in pp directory
        self.get_USPP(self.workpath_pppath)
        
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
        if not pp_files:
            qe_USPP = os.path.abspath("/public/home/mayuan/POT/qe-pp/all_pbe_UPF_v1.5")
            self.final_choosed_pp = []
            for species in self.species:
                species_name = species.name.lower()
                ppfiles       = os.listdir(qe_USPP)
                targetppfiles = filter(lambda file: species_name in file.lower(), ppfiles)
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

    # not split mode
    def get_q_from_scfout(self, dir):
        if not os.path.exists(dir):
            raise FileExistsError ("scf.out didn't exist!")
        content = open(dir, "r").readlines()
        def find_k(item):
            if re.search(r"k\(\s*\d+\)\s*=\s*", item):
                return item


        self.q_coordinate_list = []
        self.q_weight_list     = []
        self.q_total_amount    = self.q1 * self.q2 * self.q3

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

    def get_q_from_dyn0(self, dir):
        if not os.path.exists(dir):
            raise FileExistsError ("dyn0 file doesn't exist!")
        content = open(dir, "r").readlines()
        _q_total_amount = content[0].strip("\n").split()
        q_total_amount  = list(map(int, _q_total_amount))
        if q_total_amount != [self.q1, self.q2, self.q3]:
            raise ValueError ("q points set wrong")
        def find_q(item):
            if re.search(r"E\+", item):
                return item
        self.q_total_amount           = self.q1 * self.q2 * self.q3
        self.q_non_irreducible_amount = content[1]
        _q_coordinate_list            = list(filter(find_q, content))
        self.q_coordinate_list        = [q_string.strip("\n").split() for q_string in _q_coordinate_list]

        return  self.q_total_amount, self.q_non_irreducible_amount, \
                self.q_coordinate_list
    
   
    # split mode1
   
        split_ph = os.path.join(many_split_ph_dirs, "split_ph.in")
        with open(split_ph, "w") as qe:
            qe.write("Electron-phonon coefficients for {}                \n".format(self.system_name))                                    
            qe.write(" &inputph                                          \n")      
            qe.write("  tr2_ph=1.0d-16,                                  \n")              
            qe.write("  prefix='{}',                                     \n".format(self.system_name))                
            qe.write("  fildvscf='{}.dv',                                \n".format(self.system_name))                     
            qe.write("  electron_phonon='interpolated',                  \n")                              
            qe.write("  el_ph_sigma=0.005,                               \n")                 
            qe.write("  el_ph_nsigma=10,                                 \n")
            for i, species_name in enumerate(self.composition.keys()):
                element      = Element(species_name)
                species_mass = str(element.atomic_mass).strip("amu")
                qe.write("  amass({})={},                                \n".format(i+1, species_mass))             
            qe.write("  outdir='./tmp',                                  \n")               
            qe.write("  fildyn='{}.dyn',                                 \n".format(self.system_name))                    
            qe.write("  trans=.true.,                                    \n")            
            qe.write("  ldisp=.false.,                                   \n")            
            qe.write("/                                                  \n")
            qe.write(" {:<30} {:<30} {:<30}                              \n".format(q3[0], q3[1], q3[2]))
    
   
        
        split_ph_in = "split_ph" + str(start_q) + "-" + str(last_q) + ".in"
        split_ph_path = os.path.join(dir, split_ph_in)
        with open(split_ph_path, "w") as qe:
            qe.write("Electron-phonon coefficients for {}                \n".format(self.system_name))                                    
            qe.write(" &inputph                                          \n")      
            qe.write("  tr2_ph=1.0d-16,                                  \n")              
            qe.write("  prefix='{}',                                     \n".format(self.system_name))                
            qe.write("  fildvscf='{}.dv',                                \n".format(self.system_name))                     
            qe.write("  electron_phonon='interpolated',                  \n")                              
            qe.write("  el_ph_sigma=0.005,                               \n")                 
            qe.write("  el_ph_nsigma=10,                                 \n")
            for i, species_name in enumerate(self.composition.keys()):
                element      = Element(species_name)
                species_mass = str(element.atomic_mass).strip("amu")
                qe.write("  amass({})={},                                \n".format(i+1, species_mass))             
            qe.write("  outdir='./tmp',                                  \n")               
            qe.write("  fildyn='{}.dyn',                                 \n".format(self.system_name))                    
            qe.write("  trans=.true.,                                    \n")            
            qe.write("  ldisp=.true.,                                    \n")
            qe.write("  nq1={},nq2={},nq3={},                            \n".format(self.q1, self.q2, self.q3))                 
            qe.write("  start_q={}                                       \n".format(start_q)) 
            qe.write("  last_q={}                                        \n".format(last_q)) 
            qe.write("/                                                  \n")

    def merge(self, dir):
        elph_dir_path = os.path.join(dir, "elph_dir")
        if not os.path.exists(elph_dir_path):
            os.makedirs(elph_dir_path)

        for i in range(int(self.q_non_irreducible_amount)):
            src_elph   = os.path.join(dir, str(i+1), elph_dir_path, "elph.inp_lambda.1"        )
            dst_elph   = os.path.join(elph_dir_path,                "elph.inp_lambda."+str(i+1))
            shutil.copy(src_elph, dst_elph)
            logger.info(f"elph.inp_lambda.1 copy finished \n {dst_elph}")

            src_dyn    = os.path.join(dir, str(i+1), self.system_name+".dyn")
            dst_dyn    = os.path.join(dir,           self.system_name+".dyn"+str(i+1))
            shutil.copy(src_dyn, dst_dyn)
            logger.info(f"{self.system_name}.dyn copy finished \n {dst_dyn}")

            for j in range(51, 61):
                src_a2Fq2r = os.path.join(dir, str(i+1), elph_dir_path, "a2Fq2r."+str(j)+".1")
                dst_a2Fq2r = os.path.join(elph_dir_path,                "a2Fq2r."+str(j)+"."+str(i+1))
                shutil.copy(src_a2Fq2r, dst_a2Fq2r)
                logger.info(f"a2Fq2r.{str(j)}.1 copy finished \n {dst_dyn}")

        
    

    