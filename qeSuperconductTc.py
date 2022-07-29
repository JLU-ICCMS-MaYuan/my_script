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

from qe_submitjob import qe_submitjob

logging.basicConfig(
    level = logging.INFO, 
    format='%(asctime)s | %(name)s | %(levelname)s ---- %(message)s'
    )
logger = logging.getLogger(__name__)

class qe_superconduct_workflow:

    def __init__(
        self,
        input_file_path, 
        work_path,
        **kwargs,
        ):
        self.work_path            = work_path
        self.input_file_path      = input_file_path
        relax_ase                 = read(self.input_file_path)
        self.struct               = AseAtomsAdaptor.get_structure(relax_ase)
        self.get_struct_info(self.struct)
        self.work_underpressure   = None
        # self.run_relax   , self.run_scfFit     ,                      \
        # self.run_scf     , self.run_ph_no_split, self.run_ph_split  , \
        # self.run_q2r     , self.run_matdyn     , self.run_matdyn_dos, \
        # self.run_lambda = False, False, False, False, False, False, False, False, False
        if kwargs:
            for key, value in kwargs.items():
                if key == "pressure":
                    self.pressure = value
                    self.work_underpressure = os.path.join(self.work_path, str(self.pressure))
                if key == "kpoints_dense":
                    self.k1_dense , self.k2_dense , self.k3_dense  = value
                    self.k1_sparse, self.k2_sparse, self.k3_sparse = [kp/2 for kp in kwargs["kpoints_dense"]]
                    self.q1,        self.q2,        self.q3        = [kp/4 for kp in kwargs["kpoints_dense"]]
                if key == "kpoints_sparse":
                    self.k1_sparse, self.k2_sparse, self.k3_sparse = value
                    self.k1_dense , self.k2_dense , self.k3_dense  = [kp*2 for kp in kwargs["kpoints_sparse"]]
                    self.q1,        self.q2,        self.q3        = [kp/2 for kp in kwargs["kpoints_sparse"]]
                if key == "qpoints":
                    self.q1, self.q2, self.q3 = value 
                    self.k1_dense , self.k2_dense , self.k3_dense  = [kp*4 for kp in kwargs["qpoints"]]
                    self.k1_sparse, self.k2_sparse, self.k3_sparse = [kp*2 for kp in kwargs["qpoints"]]
                if key == "run_mode":
                    self.run_mode = value

        if self.work_underpressure is None:
            self.work_underpressure = self.work_path

        ####################### create work_directory ###################
        if not os.path.exists(self.work_underpressure):
            os.makedirs(self.work_underpressure)
        #######################   done work_directory ###################

        ############################ prepare pp #########################
        logger.info(f"create pp dir in {self.work_underpressure}")
        if self.input_file_path is None:
            logger.error("please specify the inputfile *.vasp or POSCAR")
            raise ValueError ("please specify the inputfile *.vasp or POSCAR")
        self.workpath_pppath = os.path.abspath(os.path.join(self.work_underpressure, "pp"))
        if not os.path.exists(self.workpath_pppath):
            os.makedirs(self.workpath_pppath)
        self.get_USPP(self.workpath_pppath)
        ############################# done pp ##########################

        # write submit task scripts
        self.submit = qe_submitjob(
            submit_path=self.work_underpressure,
            submit_job_system=self.submit_job_system,
            running_mode=self.running_mode,
            system_name=self.system_name
        )
        
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
            qe_USPP = os.path.abspath("/work/home/mayuan/POT/qe-pp/all_pbe_UPF_v1.5/")
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

    def write_relax_in(self):
        relax_in = os.path.join(self.work_underpressure, "relax.in")
        with open(relax_in, "w") as qe:
            qe.write("&CONTROL\n")
            qe.write(" calculation='vc-relax',         \n")
            qe.write(" restart_mode='from_scratch',    \n")
            qe.write(" prefix='{}',                    \n".format(self.system_name))
            qe.write(" pseudo_dir='{}',                \n".format(self.workpath_pppath))
            qe.write(" verbosity = 'high',             \n")  # 
            qe.write(" outdir='./tmp',                 \n")
            qe.write(" forc_conv_thr = 1.0d-5,         \n")
            qe.write(" etot_conv_thr = 1.0d-7,         \n")
            qe.write(" wf_collect = .true.,            \n")
            qe.write(" tstress = .true.,               \n")
            qe.write(" tprnfor = .true.,               \n")
            qe.write(" nstep = 2000,                    \n")
            qe.write("/\n")

            qe.write("&SYSTEM\n")
            qe.write(" ibrav=0,                        \n")  # 设置ibrav=0，这时需要在输入文件中写入CELL_PARAMETERS，即CELL的基矢量. alat bohr angstrom alat 由 celldm(1)或A定义的晶格常数单位
            qe.write(" nat={},                         \n".format(self.all_atoms_quantity))
            qe.write(" ntyp={},                        \n".format(self.species_quantity))
            qe.write(" occupations = 'smearing',       \n")
            # qe.write(" smearing = 'methfessel-paxton'  \n")
            qe.write(" smearing = 'gauss',             \n ")
            qe.write(" degauss = 0.005                 \n")
            qe.write(" ecutwfc = 60,                   \n")
            qe.write(" ecutrho = 720,                  \n")
            qe.write("/\n")

            qe.write("&ELECTRONS\n")
            qe.write(" diagonalization = 'david'          \n")
            qe.write(" conv_thr = 1.0d-8,              \n")
            qe.write(" mixing_beta = 0.7,              \n")
            qe.write("/\n")

            qe.write("&ions                            \n")
            qe.write(" ion_dynamics = 'bfgs',          \n")
            qe.write("/\n")

            qe.write("&cell                            \n")
            qe.write(" cell_dynamics = 'bfgs',         \n")
            qe.write(" press = {},                     \n".format(self.pressure*10))
            qe.write(" press_conv_thr = 0.01,          \n")
            qe.write("/\n")

            qe.write("ATOMIC_SPECIES                   \n")
            for species_name in self.composition.keys():
                for species_pseudo in self.final_choosed_pp:
                    if species_name.lower() in species_pseudo.lower():
                        element      = Element(species_name)
                        species_mass = str(element.atomic_mass).strip("amu")
                        qe.write(" {:<5}  {:<10}  {:<50} \n".format(species_name, species_mass, species_pseudo))
            qe.write("CELL_PARAMETERS {angstrom}        \n")  # 如果选择angstrom单未，原子坐标选择分数坐标，即，ATOMIC_POSITIONS (crystal), 且不设置celldm(1). 这时alat和celldm(1)设置成v1的长度
            for cell_p in self.cell_parameters:
                cell_p = list(map(str, cell_p))
                qe.write("{:>25} {:>25} {:>25} \n".format(cell_p[0], cell_p[1], cell_p[2]))
            qe.write("ATOMIC_POSITIONS (crystal)       \n")
            for site in self.fractional_sites:
                coord = list(map(str, site.frac_coords))
                name  = re.search(r"[A-Za-z]+", str(site.species)).group()
                # 左对齐5个字符，左对齐30个字符
                qe.write("{:<5} {:<30} {:<30} {:<30} \n".format(name, coord[0], coord[1], coord[2]))
            qe.write("K_POINTS {automatic}             \n")
            qe.write(" {} {} {} 0 0 0                  \n".format(self.k1_dense , self.k2_dense , self.k3_dense))

    def write_scf_fit_in(self, dir):
        scf_fit_in = os.path.join(dir, "scf.fit.in")
        with open(scf_fit_in, "w") as qe:
            qe.write("&CONTROL\n")
            qe.write(" calculation='scf',              \n")
            qe.write(" restart_mode='from_scratch',    \n")
            qe.write(" prefix='{}',                    \n".format(self.system_name))
            qe.write(" pseudo_dir='{}',                \n".format(self.workpath_pppath))
            qe.write(" outdir='./tmp',                 \n")
            qe.write(" forc_conv_thr = 1.0d-3,         \n")
            qe.write(" etot_conv_thr = 1.0d-4,         \n")
            qe.write("/\n")

            qe.write("&SYSTEM\n")
            qe.write(" ibrav=0,                        \n")  # 设置ibrav=0，这时需要在输入文件中写入CELL_PARAMETERS，即CELL的基矢量. alat bohr angstrom alat 由 celldm(1)或A定义的晶格常数单位
            qe.write(" nat={},                         \n".format(self.all_atoms_quantity))
            qe.write(" ntyp={},                        \n".format(self.species_quantity))
            qe.write(" occupations = 'smearing',       \n")
            qe.write(" smearing = 'methfessel-paxton'  \n")
            qe.write(" degauss = 0.02                  \n")
            qe.write(" ecutwfc = 60,                   \n")
            qe.write(" ecutrho = 720,                  \n")
            qe.write(" la2F = .true.,                  \n")
            qe.write("/\n")

            qe.write("&ELECTRONS\n")
            qe.write(" conv_thr = 1.0d-9,              \n")
            qe.write(" mixing_beta = 0.8d0,            \n")
            qe.write("/\n")

            qe.write("ATOMIC_SPECIES                   \n")
            for species_name in self.composition.keys():
                for species_pseudo in self.final_choosed_pp:
                    if species_name.lower() in species_pseudo.lower():
                        element      = Element(species_name)
                        species_mass = str(element.atomic_mass).strip("amu")
                        qe.write(" {:<5}  {:<10}  {:<50} \n".format(species_name, species_mass, species_pseudo))
            qe.write("CELL_PARAMETERS {angstrom}        \n")  # 如果选择angstrom单未，原子坐标选择分数坐标，即，ATOMIC_POSITIONS (crystal), 且不设置celldm(1). 这时alat和celldm(1)设置成v1的长度
            for cell_p in self.cell_parameters:
                cell_p = list(map(str, cell_p))
                qe.write("{:>25} {:>25} {:>25} \n".format(cell_p[0], cell_p[1], cell_p[2]))
            qe.write("ATOMIC_POSITIONS (crystal)       \n")
            for site in self.fractional_sites:
                coord = list(map(str, site.frac_coords))
                name  = re.search(r"[A-Za-z]+", str(site.species)).group()
                # 左对齐5个字符，左对齐30个字符
                qe.write("{:<5} {:<30} {:<30} {:<30} \n".format(name, coord[0], coord[1], coord[2]))
            qe.write("K_POINTS {automatic}             \n")
            qe.write(" {} {} {} 0 0 0                  \n".format(self.k1_dense, self.k2_dense, self.k3_dense))

    def write_scf_in(self, dir):
        scf_in = os.path.join(dir, "scf.in")
        with open(scf_in, "w") as qe:
            qe.write("&CONTROL\n")
            qe.write(" calculation='scf',              \n")
            qe.write(" restart_mode='from_scratch',    \n")
            qe.write(" prefix='{}',                    \n".format(self.system_name))
            qe.write(" pseudo_dir='{}',                \n".format(self.workpath_pppath))
            qe.write(" outdir='./tmp',                 \n")
            qe.write(" forc_conv_thr = 1.0d-3,         \n")
            qe.write(" etot_conv_thr = 1.0d-4,         \n")
            qe.write("/\n")

            qe.write("&SYSTEM\n")
            qe.write(" ibrav=0,                        \n")  # 设置ibrav=0，这时需要在输入文件中写入CELL_PARAMETERS，即CELL的基矢量. alat bohr angstrom alat 由 celldm(1)或A定义的晶格常数单位
            qe.write(" nat={},                         \n".format(self.all_atoms_quantity))
            qe.write(" ntyp={},                        \n".format(self.species_quantity))
            qe.write(" occupations = 'smearing',       \n")
            qe.write(" smearing = 'methfessel-paxton'  \n")
            qe.write(" degauss = 0.02                  \n")
            qe.write(" ecutwfc = 60,                   \n")
            qe.write(" ecutrho = 720,                  \n")
            qe.write("/\n")

            qe.write("&ELECTRONS\n")
            qe.write(" conv_thr = 1.0d-9,              \n")
            qe.write(" mixing_beta = 0.8d0,            \n")
            qe.write("/\n")

            qe.write("ATOMIC_SPECIES                   \n")
            for species_name in self.composition.keys():
                for species_pseudo in self.final_choosed_pp:
                    if species_name.lower() in species_pseudo.lower():
                        element      = Element(species_name)
                        species_mass = str(element.atomic_mass).strip("amu")
                        qe.write(" {:<5}  {:<10}  {:<50} \n".format(species_name, species_mass, species_pseudo))
            qe.write("CELL_PARAMETERS {angstrom}        \n")  # 如果选择angstrom单未，原子坐标选择分数坐标，即，ATOMIC_POSITIONS (crystal), 且不设置celldm(1). 这时alat和celldm(1)设置成v1的长度
            for cell_p in self.cell_parameters:
                cell_p = list(map(str, cell_p))
                qe.write("{:>25} {:>25} {:>25} \n".format(cell_p[0], cell_p[1], cell_p[2]))
            qe.write("ATOMIC_POSITIONS (crystal)       \n")
            for site in self.fractional_sites:
                coord = list(map(str, site.frac_coords))
                name  = re.search(r"[A-Za-z]+", str(site.species)).group()
                # 左对齐5个字符，左对齐30个字符
                qe.write("{:<5} {:<30} {:<30} {:<30} \n".format(name, coord[0], coord[1], coord[2]))
            qe.write("K_POINTS {automatic}             \n")
            qe.write(" {} {} {} 0 0 0                  \n".format(self.k1_sparse, self.k2_sparse, self.k3_sparse))  

    # not split mode
    def write_ph_no_split_in(self, dir):
        ph_in = os.path.join(dir, "ph_no_split.in")
        with open(ph_in, "w") as qe:
            qe.write("Electron-phonon coefficients for {}                \n".format(self.system_name))                                    
            qe.write(" &inputph                                          \n")      
            qe.write("  tr2_ph=1.0d-16,                                  \n")              
            qe.write("  prefix='{}',                                     \n".format(self.system_name))                
            qe.write("  fildvscf='{}.dv',                                \n".format(self.system_name))                     
            qe.write("  electron_phonon='interpolated',                  \n")                              
            qe.write("  el_ph_sigma=0.005,                               \n")                 
            qe.write("  el_ph_nsigma=10,                                 \n")
            qe.write("  alpha_mix(1)=0.5,                                \n")  # 可以修改的更小一些, 如果用vasp计算声子谱稳定, 可以修改为0.3
            for i, species_name in enumerate(self.composition.keys()):
                element      = Element(species_name)
                species_mass = str(element.atomic_mass).strip("amu")
                qe.write("  amass({})={},                                \n".format(i+1, species_mass))             
            qe.write("  outdir='./tmp',                                  \n")               
            qe.write("  fildyn='{}.dyn',                                 \n".format(self.system_name))                    
            qe.write("  trans=.true.,                                    \n")            
            qe.write("  ldisp=.true.,                                    \n")            
            qe.write("  nq1={},nq2={},nq3={},                            \n".format(self.q1, self.q2, self.q3 ))                 
            qe.write("/                                                  \n")
    
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
    
    def write_dyn0(self, dir):
        dyn0_path = os.path.join(dir, self.system_name+".dyn0")
        with open(dyn0_path, "w") as qe:
            qe.write("{:<5} {:<5} {:<5}              \n".format(str(self.q1), str(self.q2), str(self.q3)))
            qe.write("{}                             \n".format(str(len(self.q_list))))
            for q in self.q_list:
                qe.write("{:<30}  {:<30}  {:<30}     \n".format(q[0], q[1], q[2]))

    # split mode1
    def write_split_ph_in_from_dyn0(self, many_split_ph_dirs, q3):
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
    
    # split mode2
    def write_split_ph_in_set_startlast_q(self, dir, start_q, last_q):
        
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

    def write_q2r_in(self):
        q2r_in = os.path.join(self.work_underpressure, "q2r.in")
        with open(q2r_in, "w") as qe:
            qe.write("&input                      \n")             
            qe.write("  la2F = .true.,            \n")                       
            qe.write("  zasr = 'simple',          \n")                         
            qe.write("  fildyn = '{}.dyn'         \n".format(self.system_name))                               
            qe.write("  flfrc = '{}.fc',          \n".format(self.system_name))                              
            qe.write("/                           \n")        

    def write_matdyn_in(self):
        matdyn_in = os.path.join(self.work_underpressure, "matdyn.in")
        with open(matdyn_in, "w") as qe:
            qe.write("&input                                             \n")               
            qe.write(" asr = 'simple',                                   \n")                        
            for i, species_name in enumerate(self.composition.keys()):
                element      = Element(species_name)
                species_mass = str(element.atomic_mass).strip("amu")
                qe.write(" amass({})={},                                 \n".format(i+1, species_mass))                                   
            qe.write("  flfrc = '{}.fc',                                 \n".format(self.system_name))                              
            qe.write("  flfrq='{}.freq',                                 \n".format(self.system_name))                              
            qe.write("  la2F=.true.,                                     \n")                     
            qe.write("  dos=.flase.,                                     \n")                     
            qe.write("/                                                  \n")          
            qe.write("221                                                \n")            
            qe.write("   0.000000   0.000000   0.000000   0.0            \n")

    def write_matdyn_dos_in(self):
        matdyn_dos_in = os.path.join(self.work_underpressure, "matdyn.dos.in") 
        with open(matdyn_dos_in, "w") as qe:
            qe.write("&input                                             \n")
            qe.write("   asr = 'simple',                                 \n")                                 
            for i, species_name in enumerate(self.composition.keys()):
                element      = Element(species_name)
                species_mass = str(element.atomic_mass).strip("amu")
                qe.write(" amass({})={},                                 \n".format(i+1, species_mass))                                    
            qe.write("   flfrc = '{}.fc',                                \n".format(self.system_name))                                                                                    
            qe.write("   flfrq = '{}.freq',                              \n".format(self.system_name))                                         
            qe.write("   la2F = .true.,                                  \n")                                
            qe.write("   dos = .true.,                                   \n")                               
            qe.write("   fldos = 'phonon.dos',                           \n")                                       
            qe.write("   nk1=8, nk2=8, nk3=8,                            \n")                                      
            qe.write("   ndos=500,                                       \n")                           
            qe.write("/                                                  \n")                                                                

    def write_lambda_in(self):
        lambda_in      = os.path.join(self.work_underpressure, "lambda.in")
        qlist_dat_path = os.path.join(self.work_underpressure, "qlist.dat")
        nqs_dat_path   = os.path.join(self.work_underpressure, "nqs.dat")
        elph_dir_path  = os.path.join(self.work_underpressure, "elph_dir")
        if not os.path.exists(qlist_dat_path):
            cwd = os.getcwd()
            os.chdir(self.work_underpressure)
            os.system("sed '1,2d' *.dyn0 > qlist.dat")
            os.chdir(cwd)
        if not os.path.exists(nqs_dat_path):
            cwd = os.getcwd()
            os.chdir(self.work_underpressure)
            os.system("grep nqs q2r.out > nqs.dat")
            os.chdir(cwd)
        if not os.path.exists(elph_dir_path):
            logger.warning("There is no directory elph_dir! So the lambda.in will not be created!!!")
        else:
            a2Fq2r_elphInpLambda = os.listdir(elph_dir_path)
            elphInpLambda = sorted(list(filter(lambda x: "elph.inp_lambda" in x, a2Fq2r_elphInpLambda)))
            # prepare input data
            qlist_file = open(qlist_dat_path, "r").readlines()
            nqs_file = open(nqs_dat_path, "r").readlines()
            emax = 10;  degauss = 0.12;  smearing_method = 1
            mu = 0.01
            if len(qlist_file) == len(nqs_file) == len(elphInpLambda):
                q_number = len(nqs_file)
            else:
                logger.error("q number is wrong. The q number in qlist.dat is not matched with nqs.dat")
            with open(lambda_in, "w") as qe:
                qe.write("{:<10} {:<10} {:<10}                           \n".format(str(emax), str(degauss), str(smearing_method)))
                qe.write("{:<10}                                         \n".format(str(q_number)))
                for qlist, nq in zip(qlist_file, nqs_file):
                    qlist = qlist.strip("\n")
                    nq    =    nq.strip("\n").split("=")[-1]
                    qe.write(" {}  {}                                    \n".format(str(qlist), str(nq)))
                for elph in elphInpLambda:
                    qe.write(" {}                                        \n".format(os.path.join("elph_dir", elph)))
                qe.write("{}                                             \n".format(str(mu)))


        
    

    