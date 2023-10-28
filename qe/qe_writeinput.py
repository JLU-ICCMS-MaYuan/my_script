import os
import re
import sys
import time
import logging
from pathlib import Path

from qe.qe_inputpara import qe_inputpara
from qe.qe_base import get_pps_for_a_element

from pymatgen.core.periodic_table import Element

logger = logging.getLogger("qe_writeinput")

class qe_writeinput:
    
    def __init__(
        self, 
        work_path: Path,
        workpath_pppath: Path,
        press: int,
        qe_workflow,
        mode: str,
        **kwargs
        ):

        self.work_path = work_path
        self.workpath_pppath = Path(workpath_pppath)
        self.press = press
        self.qe_workflow = qe_workflow
        self.mode = mode
        for key, value in kwargs.items():
            setattr(self, key, value)
        

    @classmethod
    def init_from_relaxinput(cls, other_class: qe_inputpara):
        
        self = cls(
            work_path=other_class.work_path,
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
            nbnd=other_class.nbnd,
        )
        return self

    @classmethod
    def init_from_scfinput(cls, other_class: qe_inputpara):
        
        self = cls(
            work_path=other_class.work_path,
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
            nbnd=other_class.nbnd,
        )
        return self

    @classmethod
    def init_from_phonoinput(cls, other_class: qe_inputpara):
        
        self = cls(
            work_path=other_class.work_path,
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

            # phonodos
            ndos=other_class.ndos,
        )
        return self

    @classmethod
    def init_from_eletroninput(cls, other_class: qe_inputpara):

        self = cls(
            work_path=other_class.work_path,
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

            forc_conv_thr=other_class.forc_conv_thr,
            etot_conv_thr=other_class.etot_conv_thr,
            occupations=other_class.occupations,
            smearing=other_class.smearing,
            degauss=other_class.degauss,
            ecutwfc=other_class.ecutwfc,
            ecutrho=other_class.ecutrho,
            lspinorb=other_class.lspinorb,
            noncolin=other_class.noncolin,
            diagonalization=other_class.diagonalization,
            conv_thr=other_class.conv_thr,
            mixing_beta=other_class.mixing_beta,
            electron_maxstep=other_class.electron_maxstep,
            nbnd=other_class.nbnd,

            # eletrondos 
            # lspinorb=other_class.lspinorb,
            DeltaE=other_class.DeltaE,
            emin=other_class.emin,
            emax=other_class.emax,
            kinserted=other_class.kinserted,
            path_name_coords=other_class.path_name_coords,
            nbnd=other_class.nbnd,
            ngauss=other_class.ngauss,
            pdosdegauss=other_class.pdosdegauss,
        )
        return self

    @classmethod
    def init_from_scinput(cls, other_class: qe_inputpara):

        self = cls(
            work_path=other_class.work_path,
            workpath_pppath=other_class.workpath_pppath,
            press=other_class.press,
            qe_workflow=other_class.qe_workflow,
            mode=other_class.mode,

            qirreduced=other_class.qirreduced,
            qirreduced_coords=other_class.qirreduced_coords,
            qweights=other_class.qweights,

            screen_constant=other_class.screen_constant,
            top_freq=other_class.top_freq,
            degauss=other_class.degauss,
            smearing_method=other_class.smearing_method,
            temperature_steps=other_class.temperature_steps,

            a2fdos=other_class.a2fdos,
            alpha2fdat=other_class.alpha2fdat,

            # basic parameter of control precision
            forc_conv_thr=other_class.forc_conv_thr,
            etot_conv_thr=other_class.etot_conv_thr,
            smearing=other_class.smearing,
            broaden=other_class.broaden,
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
            inputfilename = self.write_relax_in(self.work_path)
            return inputfilename
        if mode == "scffit":
            inputfilename = self.write_scf_fit_in(self.work_path)
            return inputfilename
        if mode == "scf":
            inputfilename = self.write_scf_in(self.work_path)
            return inputfilename
        if mode == "nscf":
            inputfilename = self.write_nscf_in(self.work_path)
            return inputfilename
        if mode == "nosplit":
            inputfilename = self.write_ph_no_split_in(self.work_path)
            return inputfilename
        if mode == "split_dyn0":
            inputfilename = []
            for i, q3 in enumerate(self.qirreduced_coords):
                split_ph_dir = self.work_path.joinpath(str(i+1))
                if not os.path.exists(split_ph_dir):
                    os.makedirs(split_ph_dir)
                scffit_name = self.write_scf_fit_in(split_ph_dir)
                scf_name    = self.write_scf_in(split_ph_dir)
                ph_name     = self.write_split_ph_dyn0(split_ph_dir, q3)
                inputfilename.append([scffit_name, scf_name, ph_name])
                print(f"finish input files in {i+1}")
            return inputfilename     
        if mode == "split_assignQ":
            inputfilename = []
            for i in range(self.qirreduced):
                split_ph_dir = self.work_path.joinpath(str(i+1))
                if not os.path.exists(split_ph_dir):
                    os.makedirs(split_ph_dir)
                scffit_name = self.write_scf_fit_in(split_ph_dir)
                scf_name    = self.write_scf_in(split_ph_dir)
                ph_name     = self.write_split_phassignQ(split_ph_dir, i+1, i+1)
                inputfilename.append([scffit_name, scf_name, ph_name])
                print(f"finish input files in {i+1}")         
            return inputfilename
        if mode == "q2r":
            inputfilename = self.write_q2r_in(self.work_path)
            return inputfilename
        if mode == "matdyn":
            inputfilename = self.write_matdyn_in(self.work_path)
            return inputfilename
        if mode == "eletdos":
            inputfilename = self.write_eletdos_in(self.work_path)
            return inputfilename
        if mode == "elepdos":
            inputfilename = self.write_elepdos_in(self.work_path)
            return inputfilename
        if mode == "eleband":
            inputfilename = self.write_eleband_in(self.work_path)
            return inputfilename
        if mode == "elebanddata":
            inputfilename = self.write_elebanddata_in(self.work_path)
            return inputfilename
        if mode == "phonodos":
            inputfilename = self.write_phonodos_in(self.work_path)
            return inputfilename

    def write_relax_in(self, work_directory:Path):
        inputfilename =  "relax.in"
        relax_in = work_directory.joinpath(inputfilename)
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
            if self.nbnd is not None:
                qe.write(" nbnd = .{}.,                \n".format(self.nbnd))
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
                else:
                    print("\nNote: --------------------")
                    print("    When write pseudo-potentional information, something wrong. Maybe methdo `get_pps_for_a_element` is problematic !")
                    sys.exit(1)
            qe.write("CELL_PARAMETERS {angstrom}        \n")  # 如果选择angstrom单未，原子坐标选择分数坐标，即，ATOMIC_POSITIONS (crystal), 且不设置celldm(1). 这时alat和celldm(1)设置成v1的长度
            for cell_p in self.cell_parameters:
                qe.write("{}\n".format(cell_p))
            qe.write("ATOMIC_POSITIONS (crystal)       \n")
            for site in self.fractional_sites:
                qe.write("{}\n".format(site))
            qe.write("K_POINTS {automatic}             \n")
            qe.write(" {} {} {} 0 0 0                  \n".format(self.kpoints_dense[0] , self.kpoints_dense[1], self.kpoints_dense[2]))
        return inputfilename

    def write_scf_fit_in(self, work_directory:Path):
        inputfilename = "scffit.in"
        scf_fit_in = work_directory.joinpath(inputfilename)
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
            if self.nbnd is not None:
                qe.write(" nbnd = .{}.,                \n".format(self.nbnd))
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
                else:
                    print("\nNote: --------------------")
                    print("    When write pseudo-potentional information, something wrong. Maybe methdo `get_pps_for_a_element` is problematic !")
                    sys.exit(1)
            qe.write("CELL_PARAMETERS {angstrom}        \n")  # 如果选择angstrom单未，原子坐标选择分数坐标，即，ATOMIC_POSITIONS (crystal), 且不设置celldm(1). 这时alat和celldm(1)设置成v1的长度
            for cell_p in self.cell_parameters:
                qe.write("{}\n".format(cell_p))
            qe.write("ATOMIC_POSITIONS (crystal)       \n")
            for site in self.fractional_sites:
                qe.write("{}\n".format(site))
            qe.write("K_POINTS {automatic}             \n")
            qe.write(" {} {} {} 0 0 0                  \n".format(self.kpoints_dense[0] , self.kpoints_dense[1], self.kpoints_dense[2]))
        return inputfilename

    def write_scf_in(self, work_directory:Path):
        inputfilename = "scf.in"
        scf_in = work_directory.joinpath(inputfilename)
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
            if self.nbnd is not None:
                qe.write(" nbnd = .{}.,                \n".format(self.nbnd))
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
                else:
                    print("\nNote: --------------------")
                    print("    When write pseudo-potentional information, something wrong. Maybe methdo `get_pps_for_a_element` is problematic !")
                    sys.exit(1)
            qe.write("CELL_PARAMETERS {angstrom}        \n")  # 如果选择angstrom单未，原子坐标选择分数坐标，即，ATOMIC_POSITIONS (crystal), 且不设置celldm(1). 这时alat和celldm(1)设置成v1的长度
            for cell_p in self.cell_parameters:
                qe.write("{}\n".format(cell_p))
            qe.write("ATOMIC_POSITIONS (crystal)       \n")
            for site in self.fractional_sites:
                qe.write("{}\n".format(site))
            qe.write("K_POINTS {automatic}             \n")
            qe.write(" {} {} {} 0 0 0                  \n".format(self.kpoints_sparse[0] , self.kpoints_sparse[1], self.kpoints_sparse[2]))
        return inputfilename

    def write_nscf_in(self, work_directory:Path):
        inputfilename = "nscf.in"
        nscf_in = work_directory.joinpath(inputfilename)
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
            if self.nbnd is not None:
                qe.write(" nbnd = .{}.,                \n".format(self.nbnd))
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
                else:
                    print("\nNote: --------------------")
                    print("    When write pseudo-potentional information, something wrong. Maybe methdo `get_pps_for_a_element` is problematic !")
                    sys.exit(1)
            qe.write("CELL_PARAMETERS {angstrom}        \n")  # 如果选择angstrom单未，原子坐标选择分数坐标，即，ATOMIC_POSITIONS (crystal), 且不设置celldm(1). 这时alat和celldm(1)设置成v1的长度
            for cell_p in self.cell_parameters:
                qe.write("{}\n".format(cell_p))
            qe.write("ATOMIC_POSITIONS (crystal)       \n")
            for site in self.fractional_sites:
                qe.write("{}\n".format(site))
            qe.write("K_POINTS {automatic}             \n")
            qe.write(" {} {} {} 0 0 0                  \n".format(self.kpoints_dense[0], self.kpoints_dense[1], self.kpoints_dense[2]))   
        return inputfilename

    # not split mode
    def write_ph_no_split_in(self, work_directory:Path):
        inputfilename = "ph_no_split.in"
        ph_in = work_directory.joinpath(inputfilename)
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
    
    def write_dyn0(self, work_directory:Path):
        inputfilename = self.system_name+".dyn0"
        dyn0_path = work_directory.joinpath(inputfilename)
        with open(dyn0_path, "w") as qe:
            qe.write("{:<5} {:<5} {:<5}              \n".format(str(self.q1), str(self.q2), str(self.q3)))
            qe.write("{}                             \n".format(str(len(self.q_list))))
            for q in self.q_list:
                qe.write("{:<30}  {:<30}  {:<30}     \n".format(q[0], q[1], q[2]))
        return inputfilename

    # split mode1
    def write_split_ph_dyn0(self, work_directory:Path, q3):
        inputfilename = "split_ph.in"
        split_ph = work_directory.joinpath(inputfilename)
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
    def write_split_phassignQ(self, work_directory:Path, start_q, last_q):
        inputfilename = "split_ph.in"
        split_ph_path = work_directory.joinpath(inputfilename)
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
            qe.write("  outdir='./tmp',                                 \n")               
            qe.write("  fildyn='{}.dyn',                                 \n".format(self.system_name))                    
            qe.write("  trans=.true.,                                    \n")            
            qe.write("  ldisp=.true.,                                    \n")
            qe.write("  nq1={},nq2={},nq3={},                            \n".format(self.qpoints[0], self.qpoints[1], self.qpoints[2]))                 
            qe.write("  start_q={}                                       \n".format(start_q)) 
            qe.write("  last_q={}                                        \n".format(last_q)) 
            qe.write("/                                                  \n")
        return inputfilename

    def write_q2r_in(self, work_directory:Path):
        inputfilename = "q2r.in"
        q2r_in = work_directory.joinpath(inputfilename)
        with open(q2r_in, "w") as qe:
            qe.write("&input                      \n")             
            qe.write("  la2F = .true.,            \n")                       
            qe.write("  zasr = 'simple',          \n")                         
            qe.write("  fildyn = '{}.dyn'         \n".format(self.system_name))                               
            qe.write("  flfrc = '{}.fc',          \n".format(self.system_name))                              
            qe.write("/                           \n")        
        return inputfilename

    def write_matdyn_in(self, work_directory:Path):
        inputfilename = "matdyn.in"
        matdyn_in = work_directory.joinpath(inputfilename)
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

    def write_eletdos_in(self, work_directory:Path):
        inputfilename = "eletdos.in"
        eledos_in = work_directory.joinpath(inputfilename) 
        with open(eledos_in, "w") as qe:
            qe.write("&dos                       \n")
            qe.write("   prefix = '{}',          \n".format(self.system_name))                                                          
            qe.write("   outdir = './tmp',       \n")                                                                                    
            qe.write("   fildos = '{}.tdos',     \n".format(self.system_name))
            qe.write("   DeltaE = {},            \n".format(self.DeltaE))                                    
            qe.write("   emin = {},              \n".format(self.emin))
            qe.write("   emax = {},              \n".format(self.emax))
            qe.write("/                          \n")                                                                
        return inputfilename

    def write_elepdos_in(self, work_directory:Path):
        inputfilename = "elepdos.in"
        eledos_in = work_directory.joinpath(inputfilename) 
        with open(eledos_in, "w") as qe:
            qe.write("&projwfc                   \n")
            qe.write("   prefix = '{}',          \n".format(self.system_name))                                                          
            qe.write("   outdir = './tmp',       \n")                                                                                    
            qe.write("   filpdos= '{}.pdos',     \n".format(self.system_name))
            qe.write("   ngauss = {}             \n".format(self.ngauss))
            qe.write("   filproj= '{}.proj',     \n".format(self.system_name))
            qe.write("   DeltaE = {},            \n".format(self.DeltaE))   
            qe.write("   degauss = {},           \n".format(self.pdosdegauss))                                 
            qe.write("   emin = {},              \n".format(self.emin))
            qe.write("   emax = {},              \n".format(self.emax))
            qe.write("/                          \n")                                                                
        return inputfilename

    def write_eleband_in(self, work_directory:Path):
        inputfilename = "eleband.in"
        eleband_in = work_directory.joinpath(inputfilename)
        inserted_kpoints_number = self.kinserted
        with open(eleband_in, "w") as qe:
            qe.write("&CONTROL\n")
            qe.write(" calculation='bands',            \n")
            qe.write(" restart_mode='from_scratch',    \n")
            qe.write(" prefix='{}',                    \n".format(self.system_name))
            qe.write(" pseudo_dir='{}',                \n".format(str(self.workpath_pppath.absolute())))
            qe.write(" verbosity = 'high',             \n")  
            qe.write(" outdir='./tmp',                 \n")
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
            qe.write(" nbnd = {},                      \n".format(self.nbnd))
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
                else:
                    print("\nNote: --------------------")
                    print("    When write pseudo-potentional information, something wrong. Maybe methdo `get_pps_for_a_element` is problematic !")
                    sys.exit(1)
            qe.write("CELL_PARAMETERS {angstrom}        \n")  # 如果选择angstrom单未，原子坐标选择分数坐标，即，ATOMIC_POSITIONS (crystal), 且不设置celldm(1). 这时alat和celldm(1)设置成v1的长度
            for cell_p in self.cell_parameters:
                qe.write("{}\n".format(cell_p))
            qe.write("ATOMIC_POSITIONS (crystal)       \n")
            for site in self.fractional_sites:
                qe.write("{}\n".format(site))
            qe.write("K_POINTS {crystal_b}             \n")
            qe.write(" {}                              \n".format(len(self.path_name_coords)))
            for name, coord in self.path_name_coords:
                qe.write(" {:<15} {:<15} {:<15}  {}  ! {:<5}                   \n".format(str(coord[0]), str(coord[1]), str(coord[2]), str(inserted_kpoints_number), name))
        return inputfilename

    def write_elebanddata_in(self, work_directory:Path):
        inputfilename = "elebanddata.in"
        eletronband_in = work_directory.joinpath(inputfilename)
        with open(eletronband_in, "w") as f:
            f.write("&BANDS\n")
            f.write(" prefix='{}',\n".format(self.system_name))
            f.write(" outdir='./tmp',\n")
            f.write(" filband='elebanddata.dat',\n")
            f.write(" lp=.true.\n")
            f.write("/\n")
        return inputfilename
    
    def write_phonodos_in(self, work_directory:Path):
        inputfilename = "phonodos.in"
        phonodos_in = work_directory.joinpath(inputfilename) 
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

    def write_lambda_in(self, work_directory:Path, screen_constant):
        inputfilename = "lambda.in"
        lambda_in      = work_directory.joinpath(inputfilename)
        elph_dir_path  = work_directory.joinpath("elph_dir")
        if not os.path.exists(elph_dir_path):
            print("There is no directory elph_dir! So the lambda.in will not be created!!!")
            sys.exit(1)
        else:
            a2Fq2r_elphInpLambda = os.listdir(elph_dir_path)
            elphInpLambda = sorted(list(filter(lambda x: "elph.inp_lambda" in x, a2Fq2r_elphInpLambda)), key=lambda y: int(y.split('.')[-1]))
            # prepare input data
            top_freq        = self.top_freq
            broaden         = self.broaden
            smearing_method = self.smearing_method
            print("\nNote: --------------------")
            print("    Again, check to see if the four values are the same")
            print("    qirreduced_coords:{}   qweights:{}   qirreduced number:{}   elphInpLambda number:{}".format(len(self.qirreduced_coords),  len(self.qweights), int(self.qirreduced), len(elphInpLambda)))
            # time.sleep(3)
            if len(self.qirreduced_coords) == len(self.qweights) == int(self.qirreduced) == len(elphInpLambda):
                q_number = self.qirreduced
                q_coords = self.qirreduced_coords
                q_weight = self.qweights
            else:
                print("q_number is wrong. The q_number in qlist.dat is not matched with nqs.dat")
                sys.exit(1)
            with open(lambda_in, "w") as qe:
                qe.write("{:<10} {:<10} {:<10}                 \n".format(str(top_freq), str(broaden), str(smearing_method)))
                qe.write("{:<10}                               \n".format(str(q_number)))
                for qcoord, nq in zip(q_coords, q_weight):
                    qe.write(" {} {} {}  {}                    \n".format(str(qcoord[0]), str(qcoord[1]), str(qcoord[2]), str(nq)))
                for elph in elphInpLambda:
                    qe.write(" {}                              \n".format(os.path.join("elph_dir", elph)))
                qe.write("{}                                   \n".format(str(screen_constant)))
        return inputfilename

    def write_eliashberg_in(self, work_directory:Path, screen_constant):

        temperature_steps = self.temperature_steps
        inputfilename = "INPUT"
        eliashberg_in = work_directory.joinpath(inputfilename)
        with open(eliashberg_in, "w") as qe:
            qe.write("{:<10} {:<10}".format(screen_constant, temperature_steps))
        return inputfilename

    def write_alpha2f_out(self, work_directory:Path, gaussid):
        alpha2f_out = work_directory.joinpath("ALPHA2F.OUT").absolute()
        alpha2F_dat = work_directory.joinpath("alpha2F.dat").absolute()
        
        # 方法1 通过计算声子态密度, 获得a2F.dos*文件, 然后选择一个合适的文件用来计算超导
        if self.a2fdos:
            a2fdos_path = work_directory.joinpath("a2F.dos"+str(gaussid))
            if not a2fdos_path.exists():
                print(f"a2F.dos{str(gaussid)} doesn't exist !")
                sys.exit(1)
            # a2fdos_path = str(a2fdos_path)
            # alpha2f_out = str(alpha2f_out)
            os.system(f"sed '1,5d' {a2fdos_path} | sed '/lambda/d' | awk '{'{print $1/2, $2}'}' > {alpha2f_out}")
        # 方法2 不计算声子的态密度, 直接处理 alpha2F.dat 来计算超导
        else:
            if not alpha2F_dat.exists():
                print(f"{alpha2F_dat.name} doesn't exist !")
                sys.exit(1)
            # alpha2F_dat = str(alpha2F_dat)
            # alpha2f_out = str(alpha2f_out)
            awk_order   = '''sed '1,1d' %s | awk '{print $1/6579.684, $%s}' > %s''' %(alpha2F_dat, str(gaussid), alpha2f_out)
            os.system(awk_order)


