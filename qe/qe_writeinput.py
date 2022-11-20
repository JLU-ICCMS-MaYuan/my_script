import os
import re
import logging
import subprocess
from pathlib import Path

from qe.qe_inputpara import qe_inputpara
from qe.qe_base import get_pps_for_a_element

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
        self.workpath_pppath = Path(workpath_pppath)
        self.press = press
        self.qe_workflow = qe_workflow
        self.mode = mode
        for key, value in kwargs.items():
            setattr(self, key, value)
        

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

            # basic parameter of control precision
            forc_conv_thr=other_class.forc_conv_thr,
            etot_conv_thr=other_class.etot_conv_thr,
            occupations=other_class.occupations,
            smearing=other_class.smearing,
            degauss=other_class.degauss,
            ecutwfc=other_class.ecutwfc,
            ecutrho=other_class.ecutrho,
            la2F=other_class.la2F,
            lspinorb=other_class.lspinorb,
            noncolin=other_class.noncolin,
            diagonalization=other_class.diagonalization,
            conv_thr=other_class.conv_thr,
            mixing_beta=other_class.mixing_beta,
            electron_maxstep=other_class.electron_maxstep,
            press_conv_thr=other_class.press_conv_thr,
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

            # basic parameter of control precision
            forc_conv_thr=other_class.forc_conv_thr,
            etot_conv_thr=other_class.etot_conv_thr,
            occupations=other_class.occupations,
            smearing=other_class.smearing,
            degauss=other_class.degauss,
            ecutwfc=other_class.ecutwfc,
            ecutrho=other_class.ecutrho,
            la2F=other_class.la2F,
            lspinorb=other_class.lspinorb,
            noncolin=other_class.noncolin,
            diagonalization=other_class.diagonalization,
            conv_thr=other_class.conv_thr,
            mixing_beta=other_class.mixing_beta,
            electron_maxstep=other_class.electron_maxstep,
            press_conv_thr=other_class.press_conv_thr,
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

            tr2_ph=other_class.tr2_ph,
            electron_phonon=other_class.electron_phonon,
            el_ph_nsigma=other_class.el_ph_nsigma,
            el_ph_sigma=other_class.el_ph_sigma,
            alpha_mix=other_class.alpha_mix,
            kpoints_dense=other_class.kpoints_dense,
            kpoints_sparse=other_class.kpoints_sparse,
            qpoints=other_class.qpoints,

            qtot=other_class.qtot,
            qirreduced=other_class.qirreduced,
            qirreduced_coords=other_class.qirreduced_coords,
            qinserted=other_class.qinserted,
            path_name_coords=other_class.path_name_coords,

            # basic parameter of control precision
            forc_conv_thr=other_class.forc_conv_thr,
            etot_conv_thr=other_class.etot_conv_thr,
            occupations=other_class.occupations,
            smearing=other_class.smearing,
            degauss=other_class.degauss,
            ecutwfc=other_class.ecutwfc,
            ecutrho=other_class.ecutrho,
            la2F=other_class.la2F,
            lspinorb=other_class.lspinorb,
            noncolin=other_class.noncolin,
            diagonalization=other_class.diagonalization,
            conv_thr=other_class.conv_thr,
            mixing_beta=other_class.mixing_beta,
            electron_maxstep=other_class.electron_maxstep,
            press_conv_thr=other_class.press_conv_thr,
        )
        return self

    @classmethod
    def init_from_dosinput(cls, other_class: qe_inputpara):

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

            qpoints=other_class.qpoints,
            ndos=other_class.ndos,
            # lspinorb=other_class.lspinorb,
            DeltaE=other_class.DeltaE,
            emin=other_class.emin,
            emax=other_class.emax,
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

            screen_constant=other_class.screen_constant,
            top_freq=other_class.top_freq,
            deguass=other_class.deguass,
            smearing_method=other_class.smearing_method,
            temperature_points=other_class.temperature_points,
            a2F_dos=other_class.a2F_dos,
            degauss_column=other_class.degauss_column,

            # basic parameter of control precision
            forc_conv_thr=other_class.forc_conv_thr,
            etot_conv_thr=other_class.etot_conv_thr,
            smearing=other_class.smearing,
            degauss=other_class.degauss,
            ecutwfc=other_class.ecutwfc,
            ecutrho=other_class.ecutrho,
            lspinorb=other_class.lspinorb,
            noncolin=other_class.noncolin,
            diagonalization=other_class.diagonalization,
            conv_thr=other_class.conv_thr,
            mixing_beta=other_class.mixing_beta,
            press_conv_thr=other_class.press_conv_thr,
        )
        return self

    # write inputfile
    def writeinput(self, mode=None):
        if mode == None:
            mode = self.mode

        if mode == "relax-vc":
            inputfilename = self.write_relax_in()
            return inputfilename
        if mode == "scffit":
            inputfilename = self.write_scf_fit_in(self.work_underpressure)
            return inputfilename
        if mode == "scf":
            inputfilename = self.write_scf_in(self.work_underpressure)
            return inputfilename
        if mode == "nscf":
            inputfilename = self.write_nscf_in(self.work_underpressure)
            return inputfilename
        if mode =="nosplit":
            inputfilename = self.write_ph_no_split_in()
            return inputfilename
        if mode =="split_dyn0":
            inputfilename = []
            for i, q3 in enumerate(self.qirreduced_coords):
                split_ph_dir = os.path.join(self.work_underpressure, str(i+1))
                if not os.path.exists(split_ph_dir):
                    os.makedirs(split_ph_dir)
                scffit_name = self.write_scf_fit_in(split_ph_dir)
                scf_name    = self.write_scf_in(split_ph_dir)
                ph_name     = self.write_split_ph_dyn0(split_ph_dir, q3)
                inputfilename.append([scffit_name, scf_name, ph_name])
                print(f"finish input files in {i+1}")
            return inputfilename     
        if mode =="split_assignQ":
            inputfilename = []
            if self.qirreduced is not None:
                for i in range(int(self.qirreduced)):
                    ph_name     = self.write_split_phassignQ(self.work_underpressure, i+1, i+1)
                    inputfilename.append(ph_name)
                print(f"finish input files {i+1}")            
            return inputfilename
        if mode =="q2r":
            inputfilename = self.write_q2r_in()
            return inputfilename
        if mode =="matdyn":
            inputfilename = self.write_matdyn_in()
            return inputfilename
        if mode =="eletdos":
            inputfilename = self.write_eletdos_in()
            return inputfilename
        if mode =="elepdos":
            inputfilename = self.write_elepdos_in()
            return inputfilename
        if mode =="phonodos":
            inputfilename = self.write_phonodos_in()
            return inputfilename
        if mode =="McAD":
            inputfilename = self.write_lambda_in()
            return inputfilename
        if mode =="eliashberg":
            inputfilename = self.write_eliashberg_in()
            self.write_alpha2f_out()
            return inputfilename

    def write_relax_in(self):
        inputfilename =  "relax.in"
        relax_in = os.path.join(self.work_underpressure, inputfilename)
        with open(relax_in, "w") as qe:
            qe.write("&CONTROL\n")
            qe.write(" calculation='vc-relax',         \n")
            qe.write(" restart_mode='from_scratch',    \n")
            qe.write(" prefix='{}',                    \n".format(self.system_name))
            qe.write(" pseudo_dir='{}',                \n".format(str(self.workpath_pppath.absolute())))
            qe.write(" verbosity = 'high',             \n")  
            qe.write(" outdir='./tmp',                 \n")
            qe.write(" forc_conv_thr = {},             \n".format(self.forc_conv_thr))
            qe.write(" etot_conv_thr = {},             \n".format(self.etot_conv_thr))
            qe.write(" wf_collect = .true.,            \n")
            qe.write(" tstress = .true.,               \n")
            qe.write(" tprnfor = .true.,               \n")
            qe.write(" nstep = 2000,                    \n")
            qe.write("/\n")

            qe.write("&SYSTEM\n")
            qe.write(" ibrav=0,                        \n")  # 设置ibrav=0，这时需要在输入文件中写入CELL_PARAMETERS，即CELL的基矢量. alat bohr angstrom alat 由 celldm(1)或A定义的晶格常数单位
            qe.write(" nat={},                         \n".format(self.all_atoms_quantity))
            qe.write(" ntyp={},                        \n".format(self.species_quantity))
            qe.write(" occupations = '{}',             \n".format(self.occupations))
            # qe.write(" smearing = 'methfessel-paxton'  \n")
            qe.write(" smearing = '{}',                \n".format(self.smearing))
            qe.write(" degauss = {},                   \n".format(self.degauss))
            qe.write(" ecutwfc = {},                   \n".format(self.ecutwfc))
            qe.write(" ecutrho = {},                   \n".format(self.ecutrho))
            qe.write(" lspinorb = .{}.,                \n".format(self.lspinorb))
            qe.write(" noncolin = .{}.,                \n".format(self.noncolin))
            qe.write("/\n")

            qe.write("&ELECTRONS\n")
            qe.write(" diagonalization = '{}'          \n".format(self.diagonalization))
            qe.write(" conv_thr = {},                  \n".format(self.conv_thr))
            qe.write(" mixing_beta = {},               \n".format(self.mixing_beta))
            qe.write(" electron_maxstep = {},          \n".format(self.electron_maxstep))
            qe.write("/\n")

            qe.write("&ions                            \n")
            qe.write(" ion_dynamics = 'bfgs',          \n")
            qe.write("/\n")

            qe.write("&cell                            \n")
            qe.write(" cell_dynamics = 'bfgs',         \n")
            qe.write(" press = {},                     \n".format(self.press*10))
            qe.write(" press_conv_thr = {},            \n".format(self.press_conv_thr))
            qe.write("/\n")

            qe.write("ATOMIC_SPECIES                   \n")
            for species_name in self.composition.keys():
                species_pseudo = get_pps_for_a_element(species_name, self.final_choosed_pp)
                if len(species_pseudo)==1:
                    print(f"write USPP for species in relax.in: {species_pseudo[0]}") 
                    element      = Element(species_name)
                    species_mass = str(element.atomic_mass).strip("amu")
                    qe.write(" {:<5}  {:<10}  {:<50} \n".format(species_name, species_mass, species_pseudo[0]))
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
        return inputfilename

    def write_scf_fit_in(self, dir):
        inputfilename = "scffit.in"
        scf_fit_in = os.path.join(dir, inputfilename)
        with open(scf_fit_in, "w") as qe:
            qe.write("&CONTROL\n")
            qe.write(" calculation='scf',              \n")
            qe.write(" restart_mode='from_scratch',    \n")
            qe.write(" prefix='{}',                    \n".format(self.system_name))
            qe.write(" pseudo_dir='{}',                \n".format(str(self.workpath_pppath.absolute())))
            qe.write(" outdir='./tmp',                 \n")
            qe.write(" forc_conv_thr = {},             \n".format(self.forc_conv_thr))
            qe.write(" etot_conv_thr = {},             \n".format(self.etot_conv_thr))
            qe.write("/\n")

            qe.write("&SYSTEM\n")
            qe.write(" ibrav=0,                        \n")  # 设置ibrav=0，这时需要在输入文件中写入CELL_PARAMETERS，即CELL的基矢量. alat bohr angstrom alat 由 celldm(1)或A定义的晶格常数单位
            qe.write(" nat={},                         \n".format(self.all_atoms_quantity))
            qe.write(" ntyp={},                        \n".format(self.species_quantity))
            qe.write(" occupations = '{}',             \n".format(self.occupations))
            qe.write(" smearing = '{}',                \n".format(self.smearing))
            qe.write(" degauss = {},                   \n".format(self.degauss))
            qe.write(" ecutwfc = {},                   \n".format(self.ecutwfc))
            qe.write(" ecutrho = {},                   \n".format(self.ecutrho))
            qe.write(" lspinorb = .{}.,                \n".format(self.lspinorb))
            qe.write(" noncolin = .{}.,                \n".format(self.noncolin))
            qe.write(" la2F = .{}.,                    \n".format(self.la2F))
            qe.write("/\n")

            qe.write("&ELECTRONS\n")
            qe.write(" conv_thr = {},                  \n".format(self.conv_thr))
            qe.write(" mixing_beta = {},               \n".format(self.mixing_beta))
            qe.write(" electron_maxstep = {},          \n".format(self.electron_maxstep))
            qe.write("/\n")

            qe.write("ATOMIC_SPECIES                   \n")
            for species_name in self.composition.keys():
                species_pseudo = get_pps_for_a_element(species_name, self.final_choosed_pp)
                if len(species_pseudo)==1:
                    print(f"write USPP for species in scffit.in: {species_pseudo[0]}") 
                    element      = Element(species_name)
                    species_mass = str(element.atomic_mass).strip("amu")
                    qe.write(" {:<5}  {:<10}  {:<50} \n".format(species_name, species_mass, species_pseudo[0]))
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
        return inputfilename

    def write_scf_in(self, dir):
        inputfilename = "scf.in"
        scf_in = os.path.join(dir, inputfilename)
        with open(scf_in, "w") as qe:
            qe.write("&CONTROL\n")
            qe.write(" calculation='scf',              \n")
            qe.write(" restart_mode='from_scratch',    \n")
            qe.write(" prefix='{}',                    \n".format(self.system_name))
            qe.write(" pseudo_dir='{}',                \n".format(str(self.workpath_pppath.absolute())))
            qe.write(" outdir='./tmp',                 \n")
            qe.write(" forc_conv_thr = {},             \n".format(self.forc_conv_thr))
            qe.write(" etot_conv_thr = {},             \n".format(self.etot_conv_thr))
            qe.write("/\n")

            qe.write("&SYSTEM\n")
            qe.write(" ibrav=0,                        \n")  # 设置ibrav=0，这时需要在输入文件中写入CELL_PARAMETERS，即CELL的基矢量. alat bohr angstrom alat 由 celldm(1)或A定义的晶格常数单位
            qe.write(" nat={},                         \n".format(self.all_atoms_quantity))
            qe.write(" ntyp={},                        \n".format(self.species_quantity))
            qe.write(" occupations = 'smearing',       \n")
            qe.write(" smearing = '{}',                \n".format(self.smearing))
            qe.write(" degauss = {},                   \n".format(self.degauss))
            qe.write(" ecutwfc = {},                   \n".format(self.ecutwfc))
            qe.write(" ecutrho = {},                   \n".format(self.ecutrho))
            qe.write(" lspinorb = .{}.,                \n".format(self.lspinorb))
            qe.write(" noncolin = .{}.,                \n".format(self.noncolin))
            qe.write("/\n")

            qe.write("&ELECTRONS\n")
            qe.write(" conv_thr = {},                  \n".format(self.conv_thr))
            qe.write(" mixing_beta = {},               \n".format(self.mixing_beta))
            qe.write(" electron_maxstep = {},          \n".format(self.electron_maxstep))
            qe.write("/\n")

            qe.write("ATOMIC_SPECIES                   \n")
            for species_name in self.composition.keys():
                species_pseudo = get_pps_for_a_element(species_name, self.final_choosed_pp)
                if len(species_pseudo)==1:
                    print(f"write USPP for species in scf.in: {species_pseudo[0]}") 
                    element      = Element(species_name)
                    species_mass = str(element.atomic_mass).strip("amu")
                    qe.write(" {:<5}  {:<10}  {:<50} \n".format(species_name, species_mass, species_pseudo[0]))
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
        return inputfilename

    def write_nscf_in(self, dir):
        inputfilename = "nscf.in"
        nscf_in = os.path.join(dir, inputfilename)
        with open(nscf_in, "w") as qe:
            qe.write("&CONTROL\n")
            qe.write(" calculation='nscf',             \n")
            qe.write(" prefix='{}',                    \n".format(self.system_name))
            qe.write(" pseudo_dir='{}',                \n".format(str(self.workpath_pppath.absolute())))
            qe.write(" outdir='./tmp',                 \n")
            qe.write(" forc_conv_thr = {},             \n".format(self.forc_conv_thr))
            qe.write(" etot_conv_thr = {},             \n".format(self.etot_conv_thr))
            qe.write(" wf_collect=.true.,              \n")
            qe.write(" tprnfor=.true.,                 \n")
            qe.write(" tstress=.true.,                 \n")
            qe.write("/\n")

            qe.write("&SYSTEM\n")
            qe.write(" ibrav=0,                        \n")  # 设置ibrav=0，这时需要在输入文件中写入CELL_PARAMETERS，即CELL的基矢量. alat bohr angstrom alat 由 celldm(1)或A定义的晶格常数单位
            qe.write(" nat={},                         \n".format(self.all_atoms_quantity))
            qe.write(" ntyp={},                        \n".format(self.species_quantity))
            qe.write(" occupations = 'tetrahedra',     \n")
            qe.write(" ecutwfc = {},                   \n".format(self.ecutwfc))
            qe.write(" ecutrho = {},                   \n".format(self.ecutrho))
            qe.write(" lspinorb = .{}.,                \n".format(self.lspinorb))
            qe.write(" noncolin = .{}.,                \n".format(self.noncolin))
            qe.write("/\n")

            qe.write("&ELECTRONS\n")
            qe.write(" conv_thr = 1.0d-12,             \n")
            qe.write(" diagonalization = '{}'          \n".format(self.diagonalization))
            qe.write(" mixing_mode = 'plain',          \n")
            qe.write(" mixing_beta = 0.8d0,            \n")
            qe.write("/\n")

            qe.write("ATOMIC_SPECIES                   \n")
            for species_name in self.composition.keys():
                species_pseudo = get_pps_for_a_element(species_name, self.final_choosed_pp)
                if len(species_pseudo)==1:
                    print(f"write USPP for species in relax.in: {species_pseudo[0]}") 
                    element      = Element(species_name)
                    species_mass = str(element.atomic_mass).strip("amu")
                    qe.write(" {:<5}  {:<10}  {:<50} \n".format(species_name, species_mass, species_pseudo[0]))
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
            qe.write(" {} {} {} 0 0 0                  \n".format(self.kpoints_dense[0], self.kpoints_dense[1], self.kpoints_dense[2]))   
        return inputfilename

    # not split mode
    def write_ph_no_split_in(self):
        inputfilename = "ph_no_split.in"
        ph_in = os.path.join(self.work_underpressure, inputfilename)
        with open(ph_in, "w") as qe:
            qe.write("Electron-phonon coefficients for {}                \n".format(self.system_name))                                    
            qe.write(" &inputph                                          \n")      
            qe.write("  tr2_ph={},                                       \n".format(str(self.tr2_ph)))              
            qe.write("  prefix='{}',                                     \n".format(self.system_name))                
            qe.write("  fildvscf='{}.dv',                                \n".format(self.system_name))                     
            qe.write("  electron_phonon='{}',                            \n".format(self.electron_phonon))                              
            qe.write("  el_ph_sigma={},                                  \n".format(str(self.el_ph_sigma)))                
            qe.write("  el_ph_nsigma={},                                 \n".format(str(self.el_ph_nsigma)))
            qe.write("  alpha_mix(1)={},                                 \n".format(str(self.alpha_mix)))  # 可以修改的更小一些, 如果用vasp计算声子谱稳定, 可以修改为0.3
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
        return inputfilename
    
    def write_dyn0(self, dir):
        inputfilename = self.system_name+".dyn0"
        dyn0_path = os.path.join(dir, inputfilename)
        with open(dyn0_path, "w") as qe:
            qe.write("{:<5} {:<5} {:<5}              \n".format(str(self.q1), str(self.q2), str(self.q3)))
            qe.write("{}                             \n".format(str(len(self.q_list))))
            for q in self.q_list:
                qe.write("{:<30}  {:<30}  {:<30}     \n".format(q[0], q[1], q[2]))
        return inputfilename

    # split mode1
    def write_split_ph_dyn0(self, many_split_ph_dirs, q3):
        inputfilename = "split_ph.in"
        split_ph = os.path.join(many_split_ph_dirs, inputfilename)
        with open(split_ph, "w") as qe:
            qe.write("Electron-phonon coefficients for {}                \n".format(self.system_name))                                    
            qe.write(" &inputph                                          \n")      
            qe.write("  tr2_ph={},                                       \n".format(str(self.tr2_ph)))              
            qe.write("  prefix='{}',                                     \n".format(self.system_name))                
            qe.write("  fildvscf='{}.dv',                                \n".format(self.system_name))                     
            qe.write("  electron_phonon='{}',                            \n".format(self.electron_phonon))                              
            qe.write("  el_ph_sigma={},                                  \n".format(str(self.el_ph_sigma)))                
            qe.write("  el_ph_nsigma={},                                 \n".format(str(self.el_ph_nsigma)))
            qe.write("  alpha_mix(1)={},                                 \n".format(str(self.alpha_mix)))  # 可以修改的更小一些, 如果用vasp计算声子谱稳定, 可以修改为0.3
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
        return inputfilename

    # split mode2
    def write_split_phassignQ(self, dir, start_q, last_q):
        inputfilename = "splitph_" + str(start_q) + "-" + str(last_q) + ".in"
        split_ph_path = os.path.join(dir, inputfilename)
        with open(split_ph_path, "w") as qe:
            qe.write("Electron-phonon coefficients for {}                \n".format(self.system_name))                                    
            qe.write(" &inputph                                          \n")      
            qe.write("  tr2_ph={},                                       \n".format(str(self.tr2_ph)))              
            qe.write("  prefix='{}',                                     \n".format(self.system_name))                
            qe.write("  fildvscf='{}.dv',                                \n".format(self.system_name))                     
            qe.write("  electron_phonon='{}',                            \n".format(self.electron_phonon))                              
            qe.write("  el_ph_sigma={},                                  \n".format(str(self.el_ph_sigma)))                
            qe.write("  el_ph_nsigma={},                                 \n".format(str(self.el_ph_nsigma)))
            qe.write("  alpha_mix(1)={},                                 \n".format(str(self.alpha_mix)))  # 可以修改的更小一些, 如果用vasp计算声子谱稳定, 可以修改为0.3
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
        return inputfilename


    def write_q2r_in(self):
        inputfilename = "q2r.in"
        q2r_in = os.path.join(self.work_underpressure, inputfilename)
        with open(q2r_in, "w") as qe:
            qe.write("&input                      \n")             
            qe.write("  la2F = .true.,            \n")                       
            qe.write("  zasr = 'simple',          \n")                         
            qe.write("  fildyn = '{}.dyn'         \n".format(self.system_name))                               
            qe.write("  flfrc = '{}.fc',          \n".format(self.system_name))                              
            qe.write("/                           \n")        
        return inputfilename

    def write_matdyn_in(self):
        inputfilename = "matdyn.in"
        matdyn_in = os.path.join(self.work_underpressure, inputfilename)
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
            qe.write("  dos=.false.,                                     \n")                     
            qe.write("  q_in_band_form=.true.,                           \n")                                 
            qe.write("  q_in_cryst_coord=.true.,                         \n") 
            qe.write("/                                                  \n")          
            qe.write("{}                                                 \n".format(special_qpoints_number))            
            for name, coord in self.path_name_coords:
                qe.write(" {:<15} {:<15} {:<15} {:<5}                   \n".format(str(coord[0]), str(coord[1]), str(coord[2]), str(inserted_qpoints_number)))
        return inputfilename

    def write_eletdos_in(self):
        inputfilename = "eletdos.in"
        eledos_in = os.path.join(self.work_underpressure, inputfilename) 
        with open(eledos_in, "w") as qe:
            qe.write("&dos                                               \n")
            qe.write("   prefix = '{}',                                  \n".format(self.system_name))                                                          
            qe.write("   outdir = './tmp',                               \n")                                                                                    
            qe.write("   fildos = '{}.tdos',                              \n".format(self.system_name))
            qe.write("   DeltaE = 0.01,                                  \n".format(self.DeltaE))                                    
            qe.write("   emin = {},                                      \n".format(self.emin))
            qe.write("   emax = {},                                      \n".format(self.emax))
            qe.write("/                                                  \n")                                                                
        return inputfilename

    def write_elepdos_in(self):
        inputfilename = "elepdos.in"
        eledos_in = os.path.join(self.work_underpressure, inputfilename) 
        with open(eledos_in, "w") as qe:
            qe.write("&projwfc                                           \n")
            qe.write("   prefix = '{}',                                  \n".format(self.system_name))                                                          
            qe.write("   outdir = './tmp',                               \n")                                                                                    
            qe.write("   filpdos= '{}',                                  \n".format(self.system_name))
            qe.write("   filproj= '{}',                                  \n".format(self.system_name))
            qe.write("   DeltaE = 0.01,                                  \n".format(self.DeltaE))                                    
            qe.write("   emin = {},                                      \n".format(self.emin))
            qe.write("   emax = {},                                      \n".format(self.emax))
            qe.write("/                                                  \n")                                                                
        return inputfilename

    def write_phonodos_in(self):
        inputfilename = "phonodos.in"
        phonodos_in = os.path.join(self.work_underpressure, inputfilename) 
        with open(phonodos_in, "w") as qe:
            qe.write("&input                                             \n")
            qe.write("   asr = 'simple',                                 \n")                                 
            for i, species_name in enumerate(self.composition.keys()):
                element      = Element(species_name)
                species_mass = str(element.atomic_mass).strip("amu")
                qe.write("   amass({})={},                                 \n".format(i+1, species_mass))                                    
            qe.write("   flfrc = '{}.fc',                                \n".format(self.system_name))                                                                                    
            qe.write("   flfrq = '{}.freq',                              \n".format(self.system_name))                                         
            qe.write("   la2F = .true.,                                  \n")                                
            qe.write("   dos = .true.,                                   \n")  # 计算声子态密度,dos必须设置为.true.                           
            qe.write("   fldos = '{}.dos',                               \n".format(self.system_name+"_phono"))                                       
            qe.write("   nk1={}, nk2={}, nk3={},                         \n".format(self.qpoints[0], self.qpoints[1], self.qpoints[2]))  # 计算态密度时要用更密的q点网格，这需设置nk1, nk2, nk3                                      
            qe.write("   ndos={},                                        \n".format(self.ndos))  # 态密度的能量刻度上的点的数目                       
            qe.write("/                                                  \n")                                                                
        return inputfilename

    def write_lambda_in(self):
        inputfilename = "lambda.in"
        lambda_in      = os.path.join(self.work_underpressure, inputfilename)
        elph_dir_path  = os.path.join(self.work_underpressure, "elph_dir")
        if not os.path.exists(elph_dir_path):
            raise FileExistsError("There is no directory elph_dir! So the lambda.in will not be created!!!")
        else:
            a2Fq2r_elphInpLambda = os.listdir(elph_dir_path)
            elphInpLambda = sorted(list(filter(lambda x: "elph.inp_lambda" in x, a2Fq2r_elphInpLambda)))
            # prepare input data
            top_freq        = self.top_freq
            deguass         = self.deguass
            smearing_method = self.smearing_method
            screen_constant = self.screen_constant
            if len(self.qirreduced_coords) == len(self.qweights) == int(self.qirreduced) == len(elphInpLambda):
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
        return inputfilename

    def write_eliashberg_in(self):

        screen_constant = self.screen_constant
        temperature_points = self.temperature_points
        inputfilename = "INPUT"
        eliashberg_in = os.path.join(self.work_underpressure, inputfilename)
        with open(eliashberg_in, "w") as qe:
            qe.write("{:<10} {:<10}".format(screen_constant, temperature_points))
        return inputfilename

    def write_alpha2f_out(self):
        alpha2f_out = Path(self.work_underpressure).joinpath("ALPHA2F.OUT").absolute()
        alpha2F_dat = Path(self.work_underpressure).joinpath("alpha2F.dat").absolute()
        # 方法1 不计算声子的态密度, 直接处理 alpha2F.dat 来计算超导
        if self.degauss_column is not None:
            if not alpha2F_dat.exists():
                raise FileExistsError(f"{alpha2F_dat.name} doesn't exist !")
            alpha2F_dat = str(alpha2F_dat)
            alpha2f_out = str(alpha2f_out)
            awk_order   = '''sed '1,1d' %s | awk '{print $1/6579.684, $%s}' > %s''' %(alpha2F_dat, str(self.degauss_column), alpha2f_out)
            os.system(awk_order)
        # 方法2 通过计算声子态密度, 获得a2F.dos*文件, 然后选择一个合适的文件用来计算超导
        elif self.a2F_dos is not None:
            a2F_dos_path = Path(self.work_underpressure).joinpath(self.a2F_dos)
            if not a2F_dos_path.exists():
                raise FileExistsError(f"{self.a2f_dos} doesn't exist !")
            a2F_dos_path = str(a2F_dos_path)
            alpha2f_out = str(alpha2f_out)
            os.system(f"sed '1,5d' {a2F_dos_path} | sed '/lambda/d' | awk '{'{print $1/2, $2}'}' > {alpha2f_out}")
        else:
            raise("You haven't set either degauss_column or a2F_dos* !!!")

