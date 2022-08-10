import os
import re
import logging
from pathlib import Path 

from qe_inputpara import qe_inputpara

from pymatgen.core.periodic_table import Element

logger = logging.getLogger("qe_writeinput")

class qe_writeinput:
    
    def __init__(self, qe_input_object, run_mode, **kwargs):
        if isinstance(qe_input_object, qe_inputpara):
            self._qe_inputpara = qe_input_object
        self.run_mode = run_mode
        self.q_non_irreducible_amount = None
        if kwargs:
            for key, value in kwargs.items():
                if key=="q_non_irreducible_amount":
                    self.q_non_irreducible_amount = value
        
        self.writeinput()



    def writeinput(self):
        if self.run_mode == "relax":
            self.write_relax_in()
        if self.run_mode == "scffit":
            self.write_scf_fit_in(self._qe_inputpara.work_underpressure)
        if self.run_mode == "scf":
            self.write_scf_in(self._qe_inputpara.work_underpressure)
        if self.run_mode =="ph_no_split":
            self.write_ph_no_split_in()
        if self.run_mode =="ph_split_from_dyn0":
            dyn0_names = list(Path(self._qe_inputpara.work_underpressure).glob("*.dyn0"))
            if len(dyn0_names)==1:
                dyn0_path = str(dyn0_names[0].absolute())
            else:
                raise FileExistsError ("exist many *.dyn0 files or no *.dyn0")
            _, _, q_coordinate_list, _ = self._qe_inputpara.get_q_from_dyn0(dyn0_path)
            for i, q3 in enumerate(q_coordinate_list):
                split_ph_dir = os.path.join(self._qe_inputpara.work_underpressure, str(i+1))
                if not os.path.exists(split_ph_dir):
                    os.makedirs(split_ph_dir)
                self.write_split_ph_in_from_dyn0(split_ph_dir, q3)
                self.write_scf_fit_in(split_ph_dir)
                self.write_scf_in(split_ph_dir)
                logger.info(f"finish input files in {i+1}")
        if self.run_mode =="ph_split_set_startlast_q":
            if self.q_non_irreducible_amount is not None:
                for i in range(self.q_non_irreducible_amount):
                    self.write_split_ph_in_set_startlast_q(
                        self._qe_inputpara.work_underpressure, 
                        i+1, 
                        i+1)
                logger.info(f"finish input files {i+1}")
        if self.run_mode =="q2r":
            self.write_q2r_in()
        if self.run_mode =="matdyn":
            self.write_matdyn_in()
        if self.run_mode =="matdyn_dos":
            self.write_matdyn_dos_in()
        if self.run_mode =="lambda":
            self.write_lambda_in()

    def write_relax_in(self):
        relax_in = os.path.join(self._qe_inputpara.work_underpressure, "relax.in")
        with open(relax_in, "w") as qe:
            qe.write("&CONTROL\n")
            qe.write(" calculation='vc-relax',         \n")
            qe.write(" restart_mode='from_scratch',    \n")
            qe.write(" prefix='{}',                    \n".format(self._qe_inputpara.system_name))
            qe.write(" pseudo_dir='{}',                \n".format(self._qe_inputpara.workpath_pppath))
            qe.write(" verbosity = 'high',             \n")  
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
            qe.write(" nat={},                         \n".format(self._qe_inputpara.all_atoms_quantity))
            qe.write(" ntyp={},                        \n".format(self._qe_inputpara.species_quantity))
            qe.write(" occupations = 'smearing',       \n")
            # qe.write(" smearing = 'methfessel-paxton'  \n")
            qe.write(" smearing = 'gauss',             \n ")
            qe.write(" degauss = 0.005                 \n")
            qe.write(" ecutwfc = 60,                   \n")
            qe.write(" ecutrho = 720,                  \n")
            qe.write("/\n")

            qe.write("&ELECTRONS\n")
            qe.write(" diagonalization = 'david'       \n")
            qe.write(" conv_thr = 1.0d-8,              \n")
            qe.write(" mixing_beta = 0.7,              \n")
            qe.write("/\n")

            qe.write("&ions                            \n")
            qe.write(" ion_dynamics = 'bfgs',          \n")
            qe.write("/\n")

            qe.write("&cell                            \n")
            qe.write(" cell_dynamics = 'bfgs',         \n")
            qe.write(" press = {},                     \n".format(self._qe_inputpara.pressure*10))
            qe.write(" press_conv_thr = 0.01,          \n")
            qe.write("/\n")

            qe.write("ATOMIC_SPECIES                   \n")
            for species_name in self._qe_inputpara.composition.keys():
                for species_pseudo in self._qe_inputpara.final_choosed_pp:
                    match_res = re.search("^"+species_name.lower()+"\_", species_pseudo)
                    if match_res is not None:
                        logger.info(f"write USPP for species in relax.in: {match_res.group()}") 
                        element      = Element(species_name)
                        species_mass = str(element.atomic_mass).strip("amu")
                        qe.write(" {:<5}  {:<10}  {:<50} \n".format(species_name, species_mass, species_pseudo))
            qe.write("CELL_PARAMETERS {angstrom}        \n")  # 如果选择angstrom单未，原子坐标选择分数坐标，即，ATOMIC_POSITIONS (crystal), 且不设置celldm(1). 这时alat和celldm(1)设置成v1的长度
            for cell_p in self._qe_inputpara.cell_parameters:
                cell_p = list(map(str, cell_p))
                qe.write("{:>25} {:>25} {:>25} \n".format(cell_p[0], cell_p[1], cell_p[2]))
            qe.write("ATOMIC_POSITIONS (crystal)       \n")
            for site in self._qe_inputpara.fractional_sites:
                coord = list(map(str, site.frac_coords))
                name  = re.search(r"[A-Za-z]+", str(site.species)).group()
                # 左对齐5个字符，左对齐30个字符
                qe.write("{:<5} {:<30} {:<30} {:<30} \n".format(name, coord[0], coord[1], coord[2]))
            qe.write("K_POINTS {automatic}             \n")
            qe.write(" {} {} {} 0 0 0                  \n".format(self._qe_inputpara.k1_dense , self._qe_inputpara.k2_dense , self._qe_inputpara.k3_dense))

    def write_scf_fit_in(self, dir):
        scf_fit_in = os.path.join(dir, "scf.fit.in")
        with open(scf_fit_in, "w") as qe:
            qe.write("&CONTROL\n")
            qe.write(" calculation='scf',              \n")
            qe.write(" restart_mode='from_scratch',    \n")
            qe.write(" prefix='{}',                    \n".format(self._qe_inputpara.system_name))
            qe.write(" pseudo_dir='{}',                \n".format(self._qe_inputpara.workpath_pppath))
            qe.write(" outdir='./tmp',                 \n")
            qe.write(" forc_conv_thr = 1.0d-3,         \n")
            qe.write(" etot_conv_thr = 1.0d-4,         \n")
            qe.write("/\n")

            qe.write("&SYSTEM\n")
            qe.write(" ibrav=0,                        \n")  # 设置ibrav=0，这时需要在输入文件中写入CELL_PARAMETERS，即CELL的基矢量. alat bohr angstrom alat 由 celldm(1)或A定义的晶格常数单位
            qe.write(" nat={},                         \n".format(self._qe_inputpara.all_atoms_quantity))
            qe.write(" ntyp={},                        \n".format(self._qe_inputpara.species_quantity))
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
            for species_name in self._qe_inputpara.composition.keys():
                for species_pseudo in self._qe_inputpara.final_choosed_pp:
                    match_res = re.search("^"+species_name.lower()+"\_", species_pseudo)
                    if match_res is not None:
                        logger.info(f"write USPP for species in scf.fit.in: {match_res.group()}") 
                        element      = Element(species_name)
                        species_mass = str(element.atomic_mass).strip("amu")
                        qe.write(" {:<5}  {:<10}  {:<50} \n".format(species_name, species_mass, species_pseudo))
            qe.write("CELL_PARAMETERS {angstrom}        \n")  # 如果选择angstrom单未，原子坐标选择分数坐标，即，ATOMIC_POSITIONS (crystal), 且不设置celldm(1). 这时alat和celldm(1)设置成v1的长度
            for cell_p in self._qe_inputpara.cell_parameters:
                cell_p = list(map(str, cell_p))
                qe.write("{:>25} {:>25} {:>25} \n".format(cell_p[0], cell_p[1], cell_p[2]))
            qe.write("ATOMIC_POSITIONS (crystal)       \n")
            for site in self._qe_inputpara.fractional_sites:
                coord = list(map(str, site.frac_coords))
                name  = re.search(r"[A-Za-z]+", str(site.species)).group()
                # 左对齐5个字符，左对齐30个字符
                qe.write("{:<5} {:<30} {:<30} {:<30} \n".format(name, coord[0], coord[1], coord[2]))
            qe.write("K_POINTS {automatic}             \n")
            qe.write(" {} {} {} 0 0 0                  \n".format(self._qe_inputpara.k1_dense, self._qe_inputpara.k2_dense, self._qe_inputpara.k3_dense))

    def write_scf_in(self, dir):
        scf_in = os.path.join(dir, "scf.in")
        with open(scf_in, "w") as qe:
            qe.write("&CONTROL\n")
            qe.write(" calculation='scf',              \n")
            qe.write(" restart_mode='from_scratch',    \n")
            qe.write(" prefix='{}',                    \n".format(self._qe_inputpara.system_name))
            qe.write(" pseudo_dir='{}',                \n".format(self._qe_inputpara.workpath_pppath))
            qe.write(" outdir='./tmp',                 \n")
            qe.write(" forc_conv_thr = 1.0d-3,         \n")
            qe.write(" etot_conv_thr = 1.0d-4,         \n")
            qe.write("/\n")

            qe.write("&SYSTEM\n")
            qe.write(" ibrav=0,                        \n")  # 设置ibrav=0，这时需要在输入文件中写入CELL_PARAMETERS，即CELL的基矢量. alat bohr angstrom alat 由 celldm(1)或A定义的晶格常数单位
            qe.write(" nat={},                         \n".format(self._qe_inputpara.all_atoms_quantity))
            qe.write(" ntyp={},                        \n".format(self._qe_inputpara.species_quantity))
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
            for species_name in self._qe_inputpara.composition.keys():
                for species_pseudo in self._qe_inputpara.final_choosed_pp:
                    match_res = re.search("^"+species_name.lower()+"\_", species_pseudo)
                    if match_res is not None:
                        logger.info(f"write USPP for species in scf.in: {match_res.group()}") 
                        element      = Element(species_name)
                        species_mass = str(element.atomic_mass).strip("amu")
                        qe.write(" {:<5}  {:<10}  {:<50} \n".format(species_name, species_mass, species_pseudo))
            qe.write("CELL_PARAMETERS {angstrom}        \n")  # 如果选择angstrom单未，原子坐标选择分数坐标，即，ATOMIC_POSITIONS (crystal), 且不设置celldm(1). 这时alat和celldm(1)设置成v1的长度
            for cell_p in self._qe_inputpara.cell_parameters:
                cell_p = list(map(str, cell_p))
                qe.write("{:>25} {:>25} {:>25} \n".format(cell_p[0], cell_p[1], cell_p[2]))
            qe.write("ATOMIC_POSITIONS (crystal)       \n")
            for site in self._qe_inputpara.fractional_sites:
                coord = list(map(str, site.frac_coords))
                name  = re.search(r"[A-Za-z]+", str(site.species)).group()
                # 左对齐5个字符，左对齐30个字符
                qe.write("{:<5} {:<30} {:<30} {:<30} \n".format(name, coord[0], coord[1], coord[2]))
            qe.write("K_POINTS {automatic}             \n")
            qe.write(" {} {} {} 0 0 0                  \n".format(self._qe_inputpara.k1_sparse, self._qe_inputpara.k2_sparse, self._qe_inputpara.k3_sparse))  

    # not split mode
    def write_ph_no_split_in(self):
        ph_in = os.path.join(self._qe_inputpara.work_underpressure, "ph_no_split.in")
        with open(ph_in, "w") as qe:
            qe.write("Electron-phonon coefficients for {}                \n".format(self._qe_inputpara.system_name))                                    
            qe.write(" &inputph                                          \n")      
            qe.write("  tr2_ph=1.0d-16,                                  \n")              
            qe.write("  prefix='{}',                                     \n".format(self._qe_inputpara.system_name))                
            qe.write("  fildvscf='{}.dv',                                \n".format(self._qe_inputpara.system_name))                     
            qe.write("  electron_phonon='interpolated',                  \n")                              
            qe.write("  el_ph_sigma=0.005,                               \n")                 
            qe.write("  el_ph_nsigma=10,                                 \n")
            qe.write("  alpha_mix(1)=0.5,                                \n")  # 可以修改的更小一些, 如果用vasp计算声子谱稳定, 可以修改为0.3
            for i, species_name in enumerate(self._qe_inputpara.composition.keys()):
                element      = Element(species_name)
                species_mass = str(element.atomic_mass).strip("amu")
                qe.write("  amass({})={},                                \n".format(i+1, species_mass))             
            qe.write("  outdir='./tmp',                                  \n")               
            qe.write("  fildyn='{}.dyn',                                 \n".format(self._qe_inputpara.system_name))                    
            qe.write("  trans=.true.,                                    \n")            
            qe.write("  ldisp=.true.,                                    \n")            
            qe.write("  nq1={},nq2={},nq3={},                            \n".format(self._qe_inputpara.q1, self._qe_inputpara.q2, self._qe_inputpara.q3 ))                 
            qe.write("/                                                  \n")
    
    
    def write_dyn0(self, dir):
        dyn0_path = os.path.join(dir, self._qe_inputpara.system_name+".dyn0")
        with open(dyn0_path, "w") as qe:
            qe.write("{:<5} {:<5} {:<5}              \n".format(str(self._qe_inputpara.q1), str(self._qe_inputpara.q2), str(self._qe_inputpara.q3)))
            qe.write("{}                             \n".format(str(len(self._qe_inputpara.q_list))))
            for q in self._qe_inputpara.q_list:
                qe.write("{:<30}  {:<30}  {:<30}     \n".format(q[0], q[1], q[2]))

    # split mode1
    def write_split_ph_in_from_dyn0(self, many_split_ph_dirs, q3):
        split_ph = os.path.join(many_split_ph_dirs, "split_ph.in")
        with open(split_ph, "w") as qe:
            qe.write("Electron-phonon coefficients for {}                \n".format(self._qe_inputpara.system_name))                                    
            qe.write(" &inputph                                          \n")      
            qe.write("  tr2_ph=1.0d-16,                                  \n")              
            qe.write("  prefix='{}',                                     \n".format(self._qe_inputpara.system_name))                
            qe.write("  fildvscf='{}.dv',                                \n".format(self._qe_inputpara.system_name))                     
            qe.write("  electron_phonon='interpolated',                  \n")                              
            qe.write("  el_ph_sigma=0.005,                               \n")                 
            qe.write("  el_ph_nsigma=10,                                 \n")
            for i, species_name in enumerate(self._qe_inputpara.composition.keys()):
                element      = Element(species_name)
                species_mass = str(element.atomic_mass).strip("amu")
                qe.write("  amass({})={},                                \n".format(i+1, species_mass))             
            qe.write("  outdir='./tmp',                                  \n")               
            qe.write("  fildyn='{}.dyn',                                 \n".format(self._qe_inputpara.system_name))                    
            qe.write("  trans=.true.,                                    \n")            
            qe.write("  ldisp=.false.,                                   \n")            
            qe.write("/                                                  \n")
            qe.write(" {:<30} {:<30} {:<30}                              \n".format(q3[0], q3[1], q3[2]))
    
    # split mode2
    def write_split_ph_in_set_startlast_q(self, dir, start_q, last_q):
        split_ph_in = "split_ph" + str(start_q) + "-" + str(last_q) + ".in"
        split_ph_path = os.path.join(dir, split_ph_in)
        with open(split_ph_path, "w") as qe:
            qe.write("Electron-phonon coefficients for {}                \n".format(self._qe_inputpara.system_name))                                    
            qe.write(" &inputph                                          \n")      
            qe.write("  tr2_ph=1.0d-16,                                  \n")              
            qe.write("  prefix='{}',                                     \n".format(self._qe_inputpara.system_name))                
            qe.write("  fildvscf='{}.dv',                                \n".format(self._qe_inputpara.system_name))                     
            qe.write("  electron_phonon='interpolated',                  \n")                              
            qe.write("  el_ph_sigma=0.005,                               \n")                 
            qe.write("  el_ph_nsigma=10,                                 \n")
            for i, species_name in enumerate(self._qe_inputpara.composition.keys()):
                element      = Element(species_name)
                species_mass = str(element.atomic_mass).strip("amu")
                qe.write("  amass({})={},                                \n".format(i+1, species_mass))             
            qe.write("  outdir='./tmp',                                  \n")               
            qe.write("  fildyn='{}.dyn',                                 \n".format(self._qe_inputpara.system_name))                    
            qe.write("  trans=.true.,                                    \n")            
            qe.write("  ldisp=.true.,                                    \n")
            qe.write("  nq1={},nq2={},nq3={},                            \n".format(self._qe_inputpara.q1, self._qe_inputpara.q2, self._qe_inputpara.q3))                 
            qe.write("  start_q={}                                       \n".format(start_q)) 
            qe.write("  last_q={}                                        \n".format(last_q)) 
            qe.write("/                                                  \n")

 
    def write_q2r_in(self):
        q2r_in = os.path.join(self._qe_inputpara.work_underpressure, "q2r.in")
        with open(q2r_in, "w") as qe:
            qe.write("&input                      \n")             
            qe.write("  la2F = .true.,            \n")                       
            qe.write("  zasr = 'simple',          \n")                         
            qe.write("  fildyn = '{}.dyn'         \n".format(self._qe_inputpara.system_name))                               
            qe.write("  flfrc = '{}.fc',          \n".format(self._qe_inputpara.system_name))                              
            qe.write("/                           \n")        

    def write_matdyn_in(self):
        matdyn_in = os.path.join(self._qe_inputpara.work_underpressure, "matdyn.in")
        special_qpoints_number  = len(self._qe_inputpara.path_name_coords)
        inserted_qpoints_number = self._qe_inputpara.inserted_points_num
        with open(matdyn_in, "w") as qe:
            qe.write("&input                                             \n")               
            qe.write(" asr = 'simple',                                   \n")                        
            for i, species_name in enumerate(self._qe_inputpara.composition.keys()):
                element      = Element(species_name)
                species_mass = str(element.atomic_mass).strip("amu")
                qe.write(" amass({})={},                                 \n".format(i+1, species_mass))                                   
            qe.write("  flfrc = '{}.fc',                                 \n".format(self._qe_inputpara.system_name))                              
            qe.write("  flfrq='{}.freq',                                 \n".format(self._qe_inputpara.system_name))                              
            qe.write("  la2F=.true.,                                     \n")                     
            qe.write("  dos=.flase.,                                     \n")                     
            qe.write("  q_in_band_form=.true.,                           \n")                                 
            qe.write("  q_in_cryst_coord=.true.,                         \n") 
            qe.write("/                                                  \n")          
            qe.write("{}                                                 \n".format(special_qpoints_number))            
            for name, coord in self._qe_inputpara.path_name_coords:
                qe.write(" {:<15} {:<15} {:<15} {:<5}                   \n".format(str(coord[0]), str(coord[1]), str(coord[2]), str(inserted_qpoints_number)))

    def write_matdyn_dos_in(self):
        matdyn_dos_in = os.path.join(self._qe_inputpara.work_underpressure, "matdyn.dos.in") 
        with open(matdyn_dos_in, "w") as qe:
            qe.write("&input                                             \n")
            qe.write("   asr = 'simple',                                 \n")                                 
            for i, species_name in enumerate(self._qe_inputpara.composition.keys()):
                element      = Element(species_name)
                species_mass = str(element.atomic_mass).strip("amu")
                qe.write(" amass({})={},                                 \n".format(i+1, species_mass))                                    
            qe.write("   flfrc = '{}.fc',                                \n".format(self._qe_inputpara.system_name))                                                                                    
            qe.write("   flfrq = '{}.freq',                              \n".format(self._qe_inputpara.system_name))                                         
            qe.write("   la2F = .true.,                                  \n")                                
            qe.write("   dos = .true.,                                   \n")                               
            qe.write("   fldos = 'phonon.dos',                           \n")                                       
            qe.write("   nk1=8, nk2=8, nk3=8,                            \n")                                      
            qe.write("   ndos=500,                                       \n")                           
            qe.write("/                                                  \n")                                                                

    def write_lambda_in(self):
        lambda_in      = os.path.join(self._qe_inputpara.work_underpressure, "lambda.in")
        elph_dir_path  = os.path.join(self._qe_inputpara.work_underpressure, "elph_dir")
        if not os.path.exists(elph_dir_path):
            logger.warning("There is no directory elph_dir! So the lambda.in will not be created!!!")
        else:
            a2Fq2r_elphInpLambda = os.listdir(elph_dir_path)
            elphInpLambda = sorted(list(filter(lambda x: "elph.inp_lambda" in x, a2Fq2r_elphInpLambda)))
            # prepare input data
            emax = 10;  degauss = 0.12;  smearing_method = 1
            mu = 0.01
            if len(self._qe_inputpara.q_coordinate_list)        == \
               len(self._qe_inputpara.q_weight_list)            == \
               int(self._qe_inputpara.q_non_irreducible_amount) == \
               len(elphInpLambda):
                q_number = self._qe_inputpara.q_non_irreducible_amount
                q_coords = self._qe_inputpara.q_coordinate_list
                q_weight = self._qe_inputpara.q_weight_list
            else:
                logger.error("q number is wrong. The q number in qlist.dat is not matched with nqs.dat")
            with open(lambda_in, "w") as qe:
                qe.write("{:<10} {:<10} {:<10}                 \n".format(str(emax), str(degauss), str(smearing_method)))
                qe.write("{:<10}                               \n".format(str(q_number)))
                for qcoord, nq in zip(q_coords, q_weight):
                    qe.write(" {} {} {}  {}                    \n".format(str(qcoord[0]), str(qcoord[1]), str(qcoord[2]), str(nq)))
                for elph in elphInpLambda:
                    qe.write(" {}                              \n".format(os.path.join("elph_dir", elph)))
                qe.write("{}                                   \n".format(str(mu)))


        
    

