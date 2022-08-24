import os
import re
import logging
from pathlib import Path

from qe_inputpara import qe_inputpara

from pymatgen.core.periodic_table import Element

logger = logging.getLogger("qe_writeinput")

class qe_writeinput:
    
    def __init__(
        self, 
        work_underpressure: Path,
        workpath_pppath: Path,
        press: int,
        qe_workflow,
        mode: str,
        **kwargs
        ):

        self.work_underpressure = work_underpressure
        self.workpath_pppath = workpath_pppath
        self.press = press
        self.qe_workflow = qe_workflow
        self.mode = mode
        for key, value in kwargs.items():
            setattr(self, key, value)

        self.writeinput()

    @classmethod
    def init_from_relaxinput(cls, other_class: qe_inputpara):
        
        self = cls(
            work_underpressure=other_class.work_underpressure,
            workpath_pppath=other_class.workpath_pppath,
            press=other_class.press,
            qe_workflow=other_class.qe_workflow,
            mode=other_class.mode,
            system_name=other_class.system_name,
            all_atoms_quantity=other_class.all_atoms_quantity,
            species_quantity=other_class.species_quantity,
            composition=other_class.composition,
            final_choosed_pp=other_class.final_choosed_pp,
            cell_parameters=other_class.cell_parameters,
            fractional_sites=other_class.fractional_sites,
            kpoints_dense=other_class.kpoints_dense,
        )
        return self

    @classmethod
    def init_from_scfinput(cls, other_class: qe_inputpara):
        
        self = cls(
            work_underpressure=other_class.work_underpressure,
            workpath_pppath=other_class.workpath_pppath,
            press=other_class.press,
            qe_workflow=other_class.qe_workflow,
            mode=other_class.mode,
            system_name=other_class.system_name,
            all_atoms_quantity=other_class.all_atoms_quantity,
            species_quantity=other_class.species_quantity,
            composition=other_class.composition,
            final_choosed_pp=other_class.final_choosed_pp,
            cell_parameters=other_class.cell_parameters,
            fractional_sites=other_class.fractional_sites,
            kpoints_dense=other_class.kpoints_dense,
            kpoints_sparse=other_class.kpoints_sparse,
        )
        return self

    @classmethod
    def init_from_phonoinput(cls, other_class: qe_inputpara):
        
        self = cls(
            work_underpressure=other_class.work_underpressure,
            workpath_pppath=other_class.workpath_pppath,
            press=other_class.press,
            qe_workflow=other_class.qe_workflow,
            mode=other_class.mode,
            
            system_name=other_class.system_name,
            all_atoms_quantity=other_class.all_atoms_quantity,
            species_quantity=other_class.species_quantity,
            composition=other_class.composition,
            final_choosed_pp=other_class.final_choosed_pp,
            cell_parameters=other_class.cell_parameters,
            fractional_sites=other_class.fractional_sites,
            kpoints_dense=other_class.kpoints_dense,
            kpoints_sparse=other_class.kpoints_sparse,
            qpoints=other_class.qpoints,

            qtot=other_class.qtot,
            qirreduced=other_class.qirreduced,
            qirreduced_coords=other_class.qirreduced_coords,
            qinserted=other_class.qinserted,
            path_name_coords=other_class.path_name_coords,
        )
        return self

    @classmethod
    def init_from_scinput(cls, other_class: qe_inputpara):

        self = cls(
            work_underpressure=other_class.work_underpressure,
            workpath_pppath=other_class.workpath_pppath,
            press=other_class.press,
            qe_workflow=other_class.qe_workflow,
            mode=other_class.mode,

            qirreduced=other_class.qirreduced,
            qirreduced_coords=other_class.qirreduced_coords,
            qweights=other_class.qweights,
            top_freq=other_class.top_freq,
            deguass=other_class.deguass,
            smearing_method=other_class.smearing_method,
            screen_constant=other_class.screen_constant
        )
        return self

    def writeinput(self):
        if self.mode == "relax-vc":
            self.write_relax_in()
        if self.mode == "scffit":
            self.write_scf_fit_in(self.work_underpressure)
        if self.mode == "scf":
            self.write_scf_in(self.work_underpressure)
        if self.mode =="nosplit":
            self.write_ph_no_split_in()
        if self.mode =="split_from_dyn0":
            for i, q3 in enumerate(self.qirreduced_coords):
                split_ph_dir = os.path.join(self.work_underpressure, str(i+1))
                if not os.path.exists(split_ph_dir):
                    os.makedirs(split_ph_dir)
                self.write_split_ph_in_from_dyn0(split_ph_dir, q3)
                self.write_scf_fit_in(split_ph_dir)
                self.write_scf_in(split_ph_dir)
                logger.info(f"finish input files in {i+1}")
        if self.mode =="split_specify_q":
            if self.qirreduced is not None:
                for i in range(int(self.qirreduced)):
                    self.write_split_ph_in_set_startlast_q(
                        self.work_underpressure, 
                        i+1, 
                        i+1)
                logger.info(f"finish input files {i+1}")
        if self.mode =="q2r":
            self.write_q2r_in()
        if self.mode =="matdyn":
            self.write_matdyn_in()
        if self.mode =="matdyn_dos":
            self.write_matdyn_dos_in()
        if self.mode =="McAD":
            self.write_lambda_in()

    def write_relax_in(self):
        relax_in = os.path.join(self.work_underpressure, "relax.in")
        with open(relax_in, "w") as qe:
            qe.write("&CONTROL\n")
            qe.write(" calculation='vc-relax',         \n")
            qe.write(" restart_mode='from_scratch',    \n")
            qe.write(" prefix='{}',                    \n".format(self.system_name))
            qe.write(" pseudo_dir='{}',                \n".format(self.workpath_pppath))
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
            qe.write(" diagonalization = 'david'       \n")
            qe.write(" conv_thr = 1.0d-8,              \n")
            qe.write(" mixing_beta = 0.7,              \n")
            qe.write("/\n")

            qe.write("&ions                            \n")
            qe.write(" ion_dynamics = 'bfgs',          \n")
            qe.write("/\n")

            qe.write("&cell                            \n")
            qe.write(" cell_dynamics = 'bfgs',         \n")
            qe.write(" press = {},                     \n".format(self.press*10))
            qe.write(" press_conv_thr = 0.01,          \n")
            qe.write("/\n")

            qe.write("ATOMIC_SPECIES                   \n")
            for species_name in self.composition.keys():
                for species_pseudo in self.final_choosed_pp:
                    match_res = re.search("^"+species_name.lower()+"\_", species_pseudo)
                    if match_res is not None:
                        logger.info(f"write USPP for species in relax.in: {match_res.group()}") 
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
            qe.write(" {} {} {} 0 0 0                  \n".format(self.kpoints_dense[0] , self.kpoints_dense[1], self.kpoints_dense[2]))

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
                    match_res = re.search("^"+species_name.lower()+"\_", species_pseudo)
                    if match_res is not None:
                        logger.info(f"write USPP for species in scf.fit.in: {match_res.group()}") 
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
            qe.write(" {} {} {} 0 0 0                  \n".format(self.kpoints_dense[0] , self.kpoints_dense[1], self.kpoints_dense[2]))

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
                    match_res = re.search("^"+species_name.lower()+"\_", species_pseudo)
                    if match_res is not None:
                        logger.info(f"write USPP for species in scf.in: {match_res.group()}") 
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
            qe.write(" {} {} {} 0 0 0                  \n".format(self.kpoints_sparse[0], self.kpoints_sparse[1], self.kpoints_sparse[2]))  

    # not split mode
    def write_ph_no_split_in(self):
        ph_in = os.path.join(self.work_underpressure, "ph_no_split.in")
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
            qe.write("  nq1={},nq2={},nq3={},                            \n".format(self.qpoints[0], self.qpoints[1], self.qpoints[2]))                 
            qe.write("/                                                  \n")
    
    
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
            qe.write("  nq1={},nq2={},nq3={},                            \n".format(self.qpoints[0], self.qpoints[1], self.qpoints[2]))                 
            qe.write("  start_q={}                                       \n".format(start_q)) 
            qe.write("  last_q={}                                        \n".format(last_q)) 
            qe.write("/                                                  \n")

 
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
        special_qpoints_number  = len(self.path_name_coords)
        inserted_qpoints_number = self.qinserted
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
            qe.write("  q_in_band_form=.true.,                           \n")                                 
            qe.write("  q_in_cryst_coord=.true.,                         \n") 
            qe.write("/                                                  \n")          
            qe.write("{}                                                 \n".format(special_qpoints_number))            
            for name, coord in self.path_name_coords:
                qe.write(" {:<15} {:<15} {:<15} {:<5}                   \n".format(str(coord[0]), str(coord[1]), str(coord[2]), str(inserted_qpoints_number)))

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
        elph_dir_path  = os.path.join(self.work_underpressure, "elph_dir")
        if not os.path.exists(elph_dir_path):
            logger.warning("There is no directory elph_dir! So the lambda.in will not be created!!!")
        else:
            a2Fq2r_elphInpLambda = os.listdir(elph_dir_path)
            elphInpLambda = sorted(list(filter(lambda x: "elph.inp_lambda" in x, a2Fq2r_elphInpLambda)))
            # prepare input data
            top_freq        = self.top_freq
            deguass         = self.deguass
            smearing_method = self.smearing_method
            screen_constant = self.screen_constant
            if len(self.qirreduced_coords)        == \
               len(self.qweights)            == \
               int(self.qirreduced) == \
               len(elphInpLambda):
                q_number = self.qirreduced
                q_coords = self.qirreduced_coords
                q_weight = self.qweights
            else:
                logger.error("q number is wrong. The q number in qlist.dat is not matched with nqs.dat")
            with open(lambda_in, "w") as qe:
                qe.write("{:<10} {:<10} {:<10}                 \n".format(str(top_freq), str(deguass), str(smearing_method)))
                qe.write("{:<10}                               \n".format(str(q_number)))
                for qcoord, nq in zip(q_coords, q_weight):
                    qe.write(" {} {} {}  {}                    \n".format(str(qcoord[0]), str(qcoord[1]), str(qcoord[2]), str(nq)))
                for elph in elphInpLambda:
                    qe.write(" {}                              \n".format(os.path.join("elph_dir", elph)))
                qe.write("{}                                   \n".format(str(screen_constant)))


        
    

