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
        qe_inputpara: qe_inputpara
        ):

        self.qe_inputpara = qe_inputpara

    # write inputfile
    def writeinput(self, mode=None):
        if mode == None:
            mode = self.qe_inputpara.mode

        if mode == "relax-vc":
            inputfilename = self.write_relax_in(self.qe_inputpara.work_path)
            return inputfilename
        if mode == "scffit":
            inputfilename = self.write_scf_fit_in(self.qe_inputpara.work_path)
            return inputfilename
        if mode == "scf":
            inputfilename = self.write_scf_in(self.qe_inputpara.work_path)
            return inputfilename
        if mode == "nscf":
            inputfilename = self.write_nscf_in(self.qe_inputpara.work_path)
            return inputfilename
        if mode == "twin":
            inputfilename = self.write_twin_in(self.qe_inputpara.work_path)
            return inputfilename
        if mode == "kel":
            inputfilename = self.write_kel_in(self.qe_inputpara.work_path)
            return inputfilename
        if mode == "lambda_mu_k":
            inputfilename = self.write_lambda_mu_k(self.qe_inputpara.work_path)
            return inputfilename
        if mode == "scdft_tc":
            inputfilename = self.write_scdft_tc(self.qe_inputpara.work_path)
            return inputfilename
        if mode == "deltaf":
            inputfilename = self.write_deltaf(self.qe_inputpara.work_path)
            return inputfilename
        if mode == "qpdos":
            inputfilename = self.write_qpdos(self.qe_inputpara.work_path)
            return inputfilename
        if mode == "nosplit":
            inputfilename = self.write_ph_no_split_in(self.qe_inputpara.work_path)
            return inputfilename
        if mode == "split_dyn0":
            inputfilename = []
            for i, q3 in enumerate(self.qe_inputpara.qirreduced_coords):
                split_ph_dir = self.qe_inputpara.work_path.joinpath(str(i+1))
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
            for i in range(self.qe_inputpara.qirreduced):
                split_ph_dir = self.qe_inputpara.work_path.joinpath(str(i+1))
                if not os.path.exists(split_ph_dir):
                    os.makedirs(split_ph_dir)
                scffit_name = self.write_scf_fit_in(split_ph_dir)
                scf_name    = self.write_scf_in(split_ph_dir)
                ph_name     = self.write_split_phassignQ(split_ph_dir, i+1, i+1)
                inputfilename.append([scffit_name, scf_name, ph_name])
                print(f"finish input files in {i+1}")         
            return inputfilename
        if mode == "q2r":
            inputfilename = self.write_q2r_in(self.qe_inputpara.work_path)
            return inputfilename
        if mode == "matdyn":
            inputfilename = self.write_matdyn_in(self.qe_inputpara.work_path)
            return inputfilename
        if mode == "phonobanddata":
            inputfilename = self.write_phonobanddata_in(self.qe_inputpara.work_path)
            return inputfilename
        if mode == "eletdos":
            inputfilename = self.write_eletdos_in(self.qe_inputpara.work_path)
            return inputfilename
        if mode == "elepdos":
            inputfilename = self.write_elepdos_in(self.qe_inputpara.work_path)
            return inputfilename
        if mode == "eleband":
            inputfilename = self.write_eleband_in(self.qe_inputpara.work_path)
            return inputfilename
        if mode == "elebanddata":
            inputfilename = self.write_elebanddata_in(self.qe_inputpara.work_path)
            return inputfilename
        if mode == "elebandprojdata":
            inputfilename = self.write_elebandprojdata_in(self.qe_inputpara.work_path)
            return inputfilename
        if mode == "phonodos":
            inputfilename = self.write_phonodos_in(self.qe_inputpara.work_path)
            return inputfilename
        if mode == "epw_energyband":
            inputfilename = self.write_epw_energyband_in(self.qe_inputpara.work_path)
            return inputfilename
        if mode == "epw_phono":
            inputfilename = self.write_epw_phono_in(self.qe_inputpara.work_path)
            return inputfilename
        if mode == "epw_elph":
            inputfilename = self.write_epw_elph_in(self.qe_inputpara.work_path)
            return inputfilename
        if mode == "epw_aniso_sc":
            inputfilename = self.write_epw_aniso_sc_in(self.qe_inputpara.work_path)
            return inputfilename
        
    def write_relax_in(self, work_directory:Path):
        inputfilename =  "relax.in"
        relax_in = work_directory.joinpath(inputfilename)
        with open(relax_in, "w") as qe:
            qe.write("&CONTROL\n")
            qe.write(" calculation='vc-relax',         \n")
            qe.write(" restart_mode='from_scratch',    \n")
            qe.write(" prefix='{}',                    \n".format(self.qe_inputpara.system_name))
            qe.write(" pseudo_dir='{}',                \n".format(str(self.qe_inputpara.workpath_pppath.absolute())))
            qe.write(" verbosity = 'high',             \n")  
            qe.write(" outdir='./tmp',                 \n")
            qe.write(" forc_conv_thr = {},             \n".format(self.qe_inputpara.forc_conv_thr))
            qe.write(" etot_conv_thr = {},             \n".format(self.qe_inputpara.etot_conv_thr))
            qe.write(" wf_collect = .true.,            \n")
            qe.write(" tstress = .true.,               \n")
            qe.write(" tprnfor = .true.,               \n")
            qe.write(" nstep = 2000,                   \n")
            qe.write("/\n")

            qe.write("&SYSTEM\n")
            qe.write(" ibrav=0,                        \n")  # 设置ibrav=0，这时需要在输入文件中写入CELL_PARAMETERS，即CELL的基矢量. alat bohr angstrom alat 由 celldm(1)或A定义的晶格常数单位
            qe.write(" nat={},                         \n".format(self.qe_inputpara.all_atoms_quantity))
            qe.write(" ntyp={},                        \n".format(self.qe_inputpara.species_quantity))
            if self.qe_inputpara.occupations == 'smearing':
                qe.write(" occupations = '{}',         \n".format(self.qe_inputpara.occupations)) 
                qe.write(" smearing = '{}',            \n".format(self.qe_inputpara.smearing))
                qe.write(" degauss = {},               \n".format(self.qe_inputpara.degauss))
            elif self.qe_inputpara.occupations == 'tetrahedra' or self.qe_inputpara.occupations == 'tetrahedra_opt':
                qe.write(" occupations = '{}',         \n".format(self.qe_inputpara.occupations)) 
            qe.write(" ecutwfc = {},                   \n".format(self.qe_inputpara.ecutwfc))
            qe.write(" ecutrho = {},                   \n".format(self.qe_inputpara.ecutrho))
            qe.write(" lspinorb = .{}.,                \n".format(self.qe_inputpara.lspinorb))
            qe.write(" noncolin = .{}.,                \n".format(self.qe_inputpara.noncolin))
            if self.qe_inputpara.nbnd is not None:
                qe.write(" nbnd = {},                \n".format(self.qe_inputpara.nbnd))
            qe.write("/\n")

            qe.write("&ELECTRONS\n")
            qe.write(" diagonalization = '{}'          \n".format(self.qe_inputpara.diagonalization))
            qe.write(" conv_thr = {},                  \n".format(self.qe_inputpara.conv_thr))
            qe.write(" mixing_beta = {},               \n".format(self.qe_inputpara.mixing_beta))
            qe.write(" electron_maxstep = {},          \n".format(self.qe_inputpara.electron_maxstep))
            qe.write("/\n")

            qe.write("&ions                            \n")
            qe.write(" ion_dynamics = 'bfgs',          \n")
            qe.write("/\n")

            qe.write("&cell                            \n")
            qe.write(" cell_dynamics = 'bfgs',         \n")
            qe.write(" press = {},                     \n".format(self.qe_inputpara.press*10))
            qe.write(" press_conv_thr = {},            \n".format(self.qe_inputpara.press_conv_thr))
            qe.write("/\n")

            qe.write("ATOMIC_SPECIES                   \n")
            for species_name in self.qe_inputpara.composition.keys():
                species_pseudo = get_pps_for_a_element(species_name, self.qe_inputpara.final_choosed_pp)
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
            for cell_p in self.qe_inputpara.cell_parameters:
                qe.write("{}\n".format(cell_p))
            qe.write("ATOMIC_POSITIONS (crystal)       \n")
            for site in self.qe_inputpara.fractional_sites:
                qe.write("{}\n".format(site))
            qe.write("K_POINTS {automatic}             \n")
            qe.write(" {} {} {} 0 0 0                  \n".format(self.qe_inputpara.kpoints_dense[0] , self.qe_inputpara.kpoints_dense[1], self.qe_inputpara.kpoints_dense[2]))
        return inputfilename

    def write_scf_fit_in(self, work_directory:Path):
        inputfilename = "scffit.in"
        scf_fit_in = work_directory.joinpath(inputfilename)
        with open(scf_fit_in, "w") as qe:
            qe.write("&CONTROL\n")
            qe.write(" calculation='scf',              \n")
            qe.write(" restart_mode='from_scratch',    \n")
            qe.write(" prefix='{}',                    \n".format(self.qe_inputpara.system_name))
            qe.write(" pseudo_dir='{}',                \n".format(str(self.qe_inputpara.workpath_pppath.absolute())))
            qe.write(" verbosity = 'high',             \n")  
            qe.write(" outdir='./tmp',                 \n")
            qe.write(" forc_conv_thr = {},             \n".format(self.qe_inputpara.forc_conv_thr))
            qe.write(" etot_conv_thr = {},             \n".format(self.qe_inputpara.etot_conv_thr))
            qe.write(" tstress=.true.,                 \n")
            qe.write(" tprnfor=.true.,                 \n")

            qe.write("/\n")

            qe.write("&SYSTEM\n")
            qe.write(" ibrav=0,                        \n")  # 设置ibrav=0，这时需要在输入文件中写入CELL_PARAMETERS，即CELL的基矢量. alat bohr angstrom alat 由 celldm(1)或A定义的晶格常数单位
            qe.write(" nat={},                         \n".format(self.qe_inputpara.all_atoms_quantity))
            qe.write(" ntyp={},                        \n".format(self.qe_inputpara.species_quantity))
            if self.qe_inputpara.occupations == 'smearing':
                qe.write(" occupations = '{}',         \n".format(self.qe_inputpara.occupations)) 
                qe.write(" smearing = '{}',            \n".format(self.qe_inputpara.smearing))
                qe.write(" degauss = {},               \n".format(self.qe_inputpara.degauss))
            elif self.qe_inputpara.occupations == 'tetrahedra' or self.qe_inputpara.occupations == 'tetrahedra_opt':
                qe.write(" occupations = '{}',         \n".format(self.qe_inputpara.occupations)) 
            qe.write(" ecutwfc = {},                   \n".format(self.qe_inputpara.ecutwfc))
            qe.write(" ecutrho = {},                   \n".format(self.qe_inputpara.ecutrho))
            qe.write(" lspinorb = .{}.,                \n".format(self.qe_inputpara.lspinorb))
            qe.write(" noncolin = .{}.,                \n".format(self.qe_inputpara.noncolin))
            qe.write(" la2F = {},                      \n".format(self.qe_inputpara.la2F))
            if self.qe_inputpara.celldm1 is not None:
                qe.write(" celldm(1) = {:<.16f},                 \n".format(self.qe_inputpara.celldm1))
            if self.qe_inputpara.nbnd is not None:
                qe.write(" nbnd = {},                \n".format(self.qe_inputpara.nbnd))
            qe.write("/\n")

            qe.write("&ELECTRONS\n")
            qe.write(" conv_thr = {},                  \n".format(self.qe_inputpara.conv_thr))
            qe.write(" mixing_beta = {},               \n".format(self.qe_inputpara.mixing_beta))
            qe.write(" electron_maxstep = {},          \n".format(self.qe_inputpara.electron_maxstep))
            qe.write("/\n")

            qe.write("ATOMIC_SPECIES                   \n")
            for species_name in self.qe_inputpara.composition.keys():
                species_pseudo = get_pps_for_a_element(species_name, self.qe_inputpara.final_choosed_pp)
                if len(species_pseudo)==1:
                    print(f"write USPP for species in scffit.in: {species_pseudo[0]}") 
                    element      = Element(species_name)
                    species_mass = str(element.atomic_mass).strip("amu")
                    qe.write(" {:<5}  {:<10}  {:<50} \n".format(species_name, species_mass, species_pseudo[0]))
                else:
                    print("\nNote: --------------------")
                    print("    When write pseudo-potentional information, something wrong. Maybe methdo `get_pps_for_a_element` is problematic !")
                    sys.exit(1)
            if "V3_Hessian" in self.qe_inputpara.input_file_path.name:
                qe.write("CELL_PARAMETERS {alat}          \n")  # 如果选择alat单位，原子坐标选择分数坐标，即，ATOMIC_POSITIONS (crystal), 且不设置celldm(1). 这时alat和celldm(1)设置成v1的长度
            else:
                qe.write("CELL_PARAMETERS {angstrom}      \n")  # 如果选择angstrom单位，原子坐标选择分数坐标，即，ATOMIC_POSITIONS (crystal), 且不设置celldm(1). 这时alat和celldm(1)设置成v1的长度
            for cell_p in self.qe_inputpara.cell_parameters:
                qe.write("{}\n".format(cell_p))
            if "V3_Hessian" in self.qe_inputpara.input_file_path.name:
                qe.write("ATOMIC_POSITIONS {alat}         \n")
            else:
                qe.write("ATOMIC_POSITIONS {crystal}      \n")  # 如果选择angstrom单位，原子坐标选择分数坐标，即，ATOMIC_POSITIONS (crystal), 且不设置celldm(1). 这时alat和celldm(1)设置成v1的长度
            for site in self.qe_inputpara.fractional_sites:
                qe.write("{}\n".format(site))
            qe.write("K_POINTS {automatic}             \n")
            qe.write(" {} {} {} 0 0 0                  \n".format(self.qe_inputpara.kpoints_dense[0] , self.qe_inputpara.kpoints_dense[1], self.qe_inputpara.kpoints_dense[2]))
        return inputfilename

    def write_scf_in(self, work_directory:Path):
        inputfilename = "scf.in"
        scf_in = work_directory.joinpath(inputfilename)
        with open(scf_in, "w") as qe:
            qe.write("&CONTROL\n")
            qe.write(" calculation='scf',              \n")
            qe.write(" restart_mode='from_scratch',    \n")
            qe.write(" prefix='{}',                    \n".format(self.qe_inputpara.system_name))
            qe.write(" pseudo_dir='{}',                \n".format(str(self.qe_inputpara.workpath_pppath.absolute())))
            qe.write(" verbosity = 'high',             \n")  
            qe.write(" outdir='./tmp',                 \n")
            qe.write(" forc_conv_thr = {},             \n".format(self.qe_inputpara.forc_conv_thr))
            qe.write(" etot_conv_thr = {},             \n".format(self.qe_inputpara.etot_conv_thr))
            qe.write(" tstress=.true.,                 \n")
            qe.write(" tprnfor=.true.,                 \n")
            qe.write("/\n")

            qe.write("&SYSTEM\n")
            qe.write(" ibrav=0,                        \n")  # 设置ibrav=0，这时需要在输入文件中写入CELL_PARAMETERS，即CELL的基矢量. alat bohr angstrom alat 由 celldm(1)或A定义的晶格常数单位
            qe.write(" nat={},                         \n".format(self.qe_inputpara.all_atoms_quantity))
            qe.write(" ntyp={},                        \n".format(self.qe_inputpara.species_quantity))
            if self.qe_inputpara.occupations == 'smearing':
                qe.write(" occupations = '{}',         \n".format(self.qe_inputpara.occupations)) 
                qe.write(" smearing = '{}',            \n".format(self.qe_inputpara.smearing))
                qe.write(" degauss = {},               \n".format(self.qe_inputpara.degauss))
            elif self.qe_inputpara.occupations == 'tetrahedra' or self.qe_inputpara.occupations == 'tetrahedra_opt':
                qe.write(" occupations = '{}',         \n".format(self.qe_inputpara.occupations))             
            qe.write(" ecutwfc = {},                   \n".format(self.qe_inputpara.ecutwfc))
            qe.write(" ecutrho = {},                   \n".format(self.qe_inputpara.ecutrho))
            qe.write(" lspinorb = .{}.,                \n".format(self.qe_inputpara.lspinorb))
            qe.write(" noncolin = .{}.,                \n".format(self.qe_inputpara.noncolin))
            if self.qe_inputpara.celldm1 is not None:
                qe.write(" celldm(1) = {:<.16f},                 \n".format(self.qe_inputpara.celldm1))
            if self.qe_inputpara.nbnd is not None:
                qe.write(" nbnd = {},                \n".format(self.qe_inputpara.nbnd))
            qe.write("/\n")

            qe.write("&ELECTRONS\n")
            qe.write(" conv_thr = {},                  \n".format(self.qe_inputpara.conv_thr))
            qe.write(" mixing_beta = {},               \n".format(self.qe_inputpara.mixing_beta))
            qe.write(" electron_maxstep = {},          \n".format(self.qe_inputpara.electron_maxstep))
            qe.write("/\n")

            qe.write("ATOMIC_SPECIES                   \n")
            for species_name in self.qe_inputpara.composition.keys():
                species_pseudo = get_pps_for_a_element(species_name, self.qe_inputpara.final_choosed_pp)
                if len(species_pseudo)==1:
                    print(f"write USPP for species in scf.in: {species_pseudo[0]}") 
                    element      = Element(species_name)
                    species_mass = str(element.atomic_mass).strip("amu")
                    qe.write(" {:<5}  {:<10}  {:<50} \n".format(species_name, species_mass, species_pseudo[0]))
                else:
                    print("\nNote: --------------------")
                    print("    When write pseudo-potentional information, something wrong. Maybe methdo `get_pps_for_a_element` is problematic !")
                    sys.exit(1)
            if "V3_Hessian" in self.qe_inputpara.input_file_path.name:
                qe.write("CELL_PARAMETERS {alat}          \n")  # 如果选择alat单位，原子坐标选择分数坐标，即，ATOMIC_POSITIONS (crystal), 且不设置celldm(1). 这时alat和celldm(1)设置成v1的长度
            else:
                qe.write("CELL_PARAMETERS {angstrom}      \n")  # 如果选择angstrom单位，原子坐标选择分数坐标，即，ATOMIC_POSITIONS (crystal), 且不设置celldm(1). 这时alat和celldm(1)设置成v1的长度
            for cell_p in self.qe_inputpara.cell_parameters:
                qe.write("{}\n".format(cell_p))
            if "V3_Hessian" in self.qe_inputpara.input_file_path.name:
                qe.write("ATOMIC_POSITIONS {alat}         \n")
            else:
                qe.write("ATOMIC_POSITIONS {crystal}      \n")  # 如果选择angstrom单位，原子坐标选择分数坐标，即，ATOMIC_POSITIONS (crystal), 且不设置celldm(1). 这时alat和celldm(1)设置成v1的长度
            for site in self.qe_inputpara.fractional_sites:
                qe.write("{}\n".format(site))
            qe.write("K_POINTS {automatic}             \n")
            qe.write(" {} {} {} 0 0 0                  \n".format(self.qe_inputpara.kpoints_sparse[0] , self.qe_inputpara.kpoints_sparse[1], self.qe_inputpara.kpoints_sparse[2]))
        return inputfilename

    def write_nscf_in(self, work_directory:Path):
        inputfilename = "nscf.in"
        nscf_in = work_directory.joinpath(inputfilename)
        with open(nscf_in, "w") as qe:
            qe.write("&CONTROL\n")
            qe.write(" calculation='nscf',             \n")
            qe.write(" prefix='{}',                    \n".format(self.qe_inputpara.system_name))
            qe.write(" pseudo_dir='{}',                \n".format(str(self.qe_inputpara.workpath_pppath.absolute())))
            qe.write(" verbosity = 'high',             \n")  
            qe.write(" outdir='./tmp',                 \n")
            qe.write(" forc_conv_thr = {},             \n".format(self.qe_inputpara.forc_conv_thr))
            qe.write(" etot_conv_thr = {},             \n".format(self.qe_inputpara.etot_conv_thr))
            qe.write(" wf_collect=.true.,              \n")
            qe.write(" tprnfor=.true.,                 \n")
            qe.write(" tstress=.true.,                 \n")
            qe.write("/\n")

            qe.write("&SYSTEM\n")
            qe.write(" ibrav=0,                        \n")  # 设置ibrav=0，这时需要在输入文件中写入CELL_PARAMETERS，即CELL的基矢量. alat bohr angstrom alat 由 celldm(1)或A定义的晶格常数单位
            qe.write(" nat={},                         \n".format(self.qe_inputpara.all_atoms_quantity))
            qe.write(" ntyp={},                        \n".format(self.qe_inputpara.species_quantity))
            if self.qe_inputpara.occupations == 'smearing':
                qe.write(" occupations = '{}',         \n".format(self.qe_inputpara.occupations)) 
                qe.write(" smearing = '{}',            \n".format(self.qe_inputpara.smearing))
                qe.write(" degauss = {},               \n".format(self.qe_inputpara.degauss))
            elif self.qe_inputpara.occupations == 'tetrahedra' or self.qe_inputpara.occupations == 'tetrahedra_opt':
                qe.write(" occupations = '{}',         \n".format(self.qe_inputpara.occupations)) 
            qe.write(" ecutwfc = {},                   \n".format(self.qe_inputpara.ecutwfc))
            qe.write(" ecutrho = {},                   \n".format(self.qe_inputpara.ecutrho))
            qe.write(" la2F = {},                      \n".format(self.qe_inputpara.la2F))
            qe.write(" lspinorb = .{}.,                \n".format(self.qe_inputpara.lspinorb))
            qe.write(" noncolin = .{}.,                \n".format(self.qe_inputpara.noncolin))
            if self.qe_inputpara.celldm1 is not None:
                qe.write(" celldm(1) = {:<.16f},                 \n".format(self.qe_inputpara.celldm1))
            if self.qe_inputpara.nbnd is not None:
                qe.write(" nbnd = {},                \n".format(self.qe_inputpara.nbnd))
            qe.write("/\n")

            qe.write("&ELECTRONS\n")
            qe.write(" conv_thr = 1.0d-12,             \n")
            qe.write(" diagonalization = '{}'          \n".format(self.qe_inputpara.diagonalization))
            qe.write(" mixing_mode = 'plain',          \n")
            qe.write(" mixing_beta = 0.8d0,            \n")
            qe.write("/\n")

            qe.write("ATOMIC_SPECIES                   \n")
            for species_name in self.qe_inputpara.composition.keys():
                species_pseudo = get_pps_for_a_element(species_name, self.qe_inputpara.final_choosed_pp)
                if len(species_pseudo)==1:
                    print(f"write USPP for species in relax.in: {species_pseudo[0]}") 
                    element      = Element(species_name)
                    species_mass = str(element.atomic_mass).strip("amu")
                    qe.write(" {:<5}  {:<10}  {:<50} \n".format(species_name, species_mass, species_pseudo[0]))
                else:
                    print("\nNote: --------------------")
                    print("    When write pseudo-potentional information, something wrong. Maybe methdo `get_pps_for_a_element` is problematic !")
                    sys.exit(1)
            if "V3_Hessian" in self.qe_inputpara.input_file_path.name:
                qe.write("CELL_PARAMETERS {alat}          \n")  # 如果选择alat单位，原子坐标选择分数坐标，即，ATOMIC_POSITIONS (crystal), 且不设置celldm(1). 这时alat和celldm(1)设置成v1的长度
            else:
                qe.write("CELL_PARAMETERS {angstrom}      \n")  # 如果选择angstrom单位，原子坐标选择分数坐标，即，ATOMIC_POSITIONS (crystal), 且不设置celldm(1). 这时alat和celldm(1)设置成v1的长度
            for cell_p in self.qe_inputpara.cell_parameters:
                qe.write("{}\n".format(cell_p))
            if "V3_Hessian" in self.qe_inputpara.input_file_path.name:
                qe.write("ATOMIC_POSITIONS {alat}         \n")
            else:
                qe.write("ATOMIC_POSITIONS {crystal}      \n")  # 如果选择angstrom单位，原子坐标选择分数坐标，即，ATOMIC_POSITIONS (crystal), 且不设置celldm(1). 这时alat和celldm(1)设置成v1的长度
            for site in self.qe_inputpara.fractional_sites:
                qe.write("{}\n".format(site))
            if self.qe_inputpara.k_automatic:
                qe.write("K_POINTS {automatic}             \n")
                qe.write(" {} {} {} 0 0 0                  \n".format(self.qe_inputpara.kpoints_dense[0], self.qe_inputpara.kpoints_dense[1], self.qe_inputpara.kpoints_dense[2]))
            else:
                qe.write("K_POINTS crystal                 \n")
                qe.write("{}\n".format(self.qe_inputpara.totpts))
                for kinfo in self.qe_inputpara.kpoints_coords:
                    qe.write(" {}    \n".format(kinfo))
        return inputfilename

    def write_twin_in(self, work_directory:Path):
        inputfilename = "twin.in"
        twin_in = work_directory.joinpath(inputfilename)
        with open(twin_in, "w") as qe:
            qe.write("&CONTROL\n")
            qe.write(" calculation='bands',             \n")
            qe.write(" prefix='{}',                    \n".format(self.qe_inputpara.system_name))
            qe.write(" pseudo_dir='{}',                \n".format(str(self.qe_inputpara.workpath_pppath.absolute())))
            qe.write(" verbosity = 'high',             \n")  
            qe.write(" outdir='./tmp',                 \n")
            qe.write(" etot_conv_thr = {},             \n".format(self.qe_inputpara.etot_conv_thr))
            qe.write("/\n")

            qe.write("&SYSTEM\n")
            qe.write(" ibrav=0,                        \n")  # 设置ibrav=0，这时需要在输入文件中写入CELL_PARAMETERS，即CELL的基矢量. alat bohr angstrom alat 由 celldm(1)或A定义的晶格常数单位
            qe.write(" nat={},                         \n".format(self.qe_inputpara.all_atoms_quantity))
            qe.write(" ntyp={},                        \n".format(self.qe_inputpara.species_quantity))
            qe.write(" ecutwfc = {},                   \n".format(self.qe_inputpara.ecutwfc))
            qe.write(" ecutrho = {},                   \n".format(self.qe_inputpara.ecutrho))
            if self.qe_inputpara.celldm1 is not None:
                qe.write(" celldm(1) = {:<.16f},                 \n".format(self.qe_inputpara.celldm1))
            if self.qe_inputpara.nbnd is not None:
                qe.write(" nbnd = {},                \n".format(self.qe_inputpara.nbnd))
            qe.write("/\n")

            qe.write("&ELECTRONS\n")
            qe.write("/\n")

            qe.write("ATOMIC_SPECIES                   \n")
            for species_name in self.qe_inputpara.composition.keys():
                species_pseudo = get_pps_for_a_element(species_name, self.qe_inputpara.final_choosed_pp)
                if len(species_pseudo)==1:
                    print(f"write USPP for species in relax.in: {species_pseudo[0]}") 
                    element      = Element(species_name)
                    species_mass = str(element.atomic_mass).strip("amu")
                    qe.write(" {:<5}  {:<10}  {:<50} \n".format(species_name, species_mass, species_pseudo[0]))
                else:
                    print("\nNote: --------------------")
                    print("    When write pseudo-potentional information, something wrong. Maybe methdo `get_pps_for_a_element` is problematic !")
                    sys.exit(1)
            if "V3_Hessian" in self.qe_inputpara.input_file_path.name:
                qe.write("CELL_PARAMETERS {alat}          \n")  # 如果选择alat单位，原子坐标选择分数坐标，即，ATOMIC_POSITIONS (crystal), 且不设置celldm(1). 这时alat和celldm(1)设置成v1的长度
            else:
                qe.write("CELL_PARAMETERS {angstrom}      \n")  # 如果选择angstrom单位，原子坐标选择分数坐标，即，ATOMIC_POSITIONS (crystal), 且不设置celldm(1). 这时alat和celldm(1)设置成v1的长度
            for cell_p in self.qe_inputpara.cell_parameters:
                qe.write("{}\n".format(cell_p))
            if "V3_Hessian" in self.qe_inputpara.input_file_path.name:
                qe.write("ATOMIC_POSITIONS {alat}         \n")
            else:
                qe.write("ATOMIC_POSITIONS {crystal}      \n")  # 如果选择angstrom单位，原子坐标选择分数坐标，即，ATOMIC_POSITIONS (crystal), 且不设置celldm(1). 这时alat和celldm(1)设置成v1的长度
            for site in self.qe_inputpara.fractional_sites:
                qe.write("{}\n".format(site))
            qe.write("K_POINTS crystal                 \n")
            qe.write("{}\n".format(self.qe_inputpara.totpts_for_Twin))
            for kinfo in self.qe_inputpara.kpoints_coords_for_Twin:
                qe.write(" {}    \n".format(kinfo))
        return inputfilename

    def write_kel_in(self, work_directory:Path):
        inputfilename = "kel.in"
        kel_in = work_directory.joinpath(inputfilename)
        with open(kel_in, "w") as qe:
            qe.write("&CONTROL\n")
            qe.write(" calculation='kel'            \n")
            qe.write(" prefix='{}',                 \n".format(self.qe_inputpara.system_name))                
            qe.write(" outdir='./tmp',              \n")               
            qe.write("/                             \n")
            qe.write("&Kel                          \n")
            qe.write(" nci=5                        \n")   
            qe.write(" laddxc=0                     \n")      
            qe.write(" lsf=1                        \n")   
            qe.write(" ecutfock=70.0                \n")
            print(len(self.qe_inputpara.qpoints))
            qe.write(" nq1={}, nq2={}, nq3={}\n".format(self.qe_inputpara.qpoints[0], self.qe_inputpara.qpoints[1], self.qe_inputpara.qpoints[2]))
            qe.write("/                             \n")
            qe.write("&SCDFT                        \n")
            qe.write(" temp = -0.1                  \n")   
            qe.write(" fbee = 1                     \n")
            qe.write(" lbee = 15                    \n") 
            qe.write(" xic = -1.0                   \n")  
            qe.write(" nmf = 10                     \n")
            qe.write(" nx = 100                     \n")
            qe.write(" ne = 50                      \n")
            qe.write(" emin = 1.0e-7                \n")     
            qe.write(" emax = 5.0                   \n")  
            qe.write(" electron_maxstep = 100       \n")
            qe.write(" conv_thr = 1.0e-15           \n")
            qe.write(" spin_fluc =.true.            \n")  
            qe.write("/                             \n")      
        return inputfilename
    
    def write_lambda_mu_k(self, work_directory:Path):
        inputfilename = "lambda_mu_k.in"
        lambda_mu_k = work_directory.joinpath(inputfilename)
        with open(lambda_mu_k, "w") as qe:
            qe.write("&CONTROL\n")
            qe.write(" calculation='lambda_mu_k'    \n")
            qe.write(" prefix='{}',                 \n".format(self.qe_inputpara.system_name))                
            qe.write(" outdir='./tmp',              \n")               
            qe.write("/                             \n")
            qe.write("&Kel                          \n")
            qe.write(" nci=5                        \n")   
            qe.write(" laddxc=0                     \n")      
            qe.write(" lsf=1                        \n")   
            qe.write(" ecutfock=70.0                \n")
            qe.write(" nq1={}, nq2={}, nq3={}\n".format(self.qe_inputpara.qpoints[0], self.qe_inputpara.qpoints[1], self.qe_inputpara.qpoints[2]))
            qe.write("/                             \n")
            qe.write("&SCDFT                        \n")
            qe.write(" temp = -0.1                  \n")   
            qe.write(" fbee = 1                     \n")
            qe.write(" lbee = 10                    \n")  # 这里不同于kel.in中的 lbee, 和scdft_tc, deltaf, qpdos中一样
            qe.write(" xic = -1.0                   \n")  
            qe.write(" nmf = 10                     \n")
            qe.write(" nx = 100                     \n")
            qe.write(" ne = 50                      \n")
            qe.write(" emin = 1.0e-7                \n")     
            qe.write(" emax = 0.7                   \n")  # 这里不同于kel.in中的 emax, 和scdft_tc, deltaf, qpdos中一样
            qe.write(" electron_maxstep = 100       \n")
            qe.write(" conv_thr = 1.0e-15           \n")
            qe.write(" spin_fluc =.true.            \n")  
            qe.write("/                             \n")      
        return inputfilename
    
    def write_scdft_tc(self, work_directory:Path):
        inputfilename = "scdft_tc.in"
        scdft_tc = work_directory.joinpath(inputfilename)
        with open(scdft_tc, "w") as qe:
            qe.write("&CONTROL\n")
            qe.write(" calculation='scdft_tc'       \n")
            qe.write(" prefix='{}',                 \n".format(self.qe_inputpara.system_name))                
            qe.write(" outdir='./tmp',              \n")               
            qe.write("/                             \n")
            qe.write("&Kel                          \n")
            qe.write(" nci=5                        \n")   
            qe.write(" laddxc=0                     \n")      
            qe.write(" lsf=1                        \n")   
            qe.write(" ecutfock=70.0                \n")
            qe.write(" nq1={}, nq2={}, nq3={}\n".format(self.qe_inputpara.qpoints[0], self.qe_inputpara.qpoints[1], self.qe_inputpara.qpoints[2]))
            qe.write("/                             \n")
            qe.write("&SCDFT                        \n")
            qe.write(" temp = -0.1                  \n")   
            qe.write(" fbee = 1                     \n")
            qe.write(" lbee = 10                    \n")  # 这里不同于kel.in中的 lbee, 和delta, lambda_mu_kf, qpdos中一样
            qe.write(" xic = -1.0                   \n")  
            qe.write(" nmf = 10                     \n")
            qe.write(" nx = 100                     \n")
            qe.write(" ne = 50                      \n")
            qe.write(" emin = 1.0e-7                \n")     
            qe.write(" emax = 0.7                   \n")  # 这里不同于kel.in中的 emax, 和delta, lambda_mu_kf, qpdos中一样
            qe.write(" electron_maxstep = 100       \n")
            qe.write(" conv_thr = 1.0e-15           \n")
            qe.write(" spin_fluc =.true.            \n")  
            qe.write("/                             \n")      
        return inputfilename
    
    def write_deltaf(self, work_directory:Path):
        inputfilename = "deltaf.in"
        deltaf = work_directory.joinpath(inputfilename)
        with open(deltaf, "w") as qe:
            qe.write("&CONTROL\n")
            qe.write(" calculation='deltaf'       \n")
            qe.write(" prefix='{}',                 \n".format(self.qe_inputpara.system_name))                
            qe.write(" outdir='./tmp',              \n")               
            qe.write("/                             \n")
            qe.write("&Kel                          \n")
            qe.write(" nci=5                        \n")   
            qe.write(" laddxc=0                     \n")      
            qe.write(" lsf=1                        \n")   
            qe.write(" ecutfock=70.0                \n")
            qe.write(" nq1={}, nq2={}, nq3={}\n".format(self.qe_inputpara.qpoints[0], self.qe_inputpara.qpoints[1], self.qe_inputpara.qpoints[2]))
            qe.write("/                             \n")
            qe.write("&SCDFT                        \n")
            qe.write(" temp = -0.1                  \n")   
            qe.write(" fbee = 1                     \n")
            qe.write(" lbee = 10                    \n")  # 这里不同于kel.in中的 lbee, 和scdft_tc, lambda_mu_kf, qpdos中一样
            qe.write(" xic = -1.0                   \n")  
            qe.write(" nmf = 10                     \n")
            qe.write(" nx = 100                     \n")
            qe.write(" ne = 50                      \n")
            qe.write(" emin = 1.0e-7                \n")     
            qe.write(" emax = 0.7                   \n")  # 这里不同于kel.in中的 emax, 和scdft_tc, lambda_mu_kf, qpdos中一样
            qe.write(" electron_maxstep = 100       \n")
            qe.write(" conv_thr = 1.0e-15           \n")
            qe.write(" spin_fluc =.true.            \n")  
            qe.write("/                             \n")      
        return inputfilename
    
    def write_qpdos(self, work_directory:Path):
        inputfilename = "qpdos.in"
        qpdos = work_directory.joinpath(inputfilename)
        with open(qpdos, "w") as qe:
            qe.write("&CONTROL\n")
            qe.write(" calculation='qpdos'          \n")
            qe.write(" prefix='{}',                 \n".format(self.qe_inputpara.system_name))                
            qe.write(" outdir='./tmp',              \n")               
            qe.write("/                             \n")
            qe.write("&Kel                          \n")
            qe.write(" nci=5                        \n")   
            qe.write(" laddxc=0                     \n")      
            qe.write(" lsf=1                        \n")   
            qe.write(" ecutfock=70.0                \n")
            qe.write(" nq1={}, nq2={}, nq3={}\n".format(self.qe_inputpara.qpoints[0], self.qe_inputpara.qpoints[1], self.qe_inputpara.qpoints[2]))
            qe.write("/                             \n")
            qe.write("&SCDFT                        \n")
            qe.write(" temp = -0.1                  \n")   
            qe.write(" fbee = 1                     \n")
            qe.write(" lbee = 10                    \n")  # 这里不同于kel.in中的 lbee, 其它都一样
            qe.write(" xic = -1.0                   \n")  
            qe.write(" nmf = 10                     \n")
            qe.write(" nx = 100                     \n")
            qe.write(" ne = 50                      \n")
            qe.write(" emin = 1.0e-7                \n")     
            qe.write(" emax = 0.7                   \n")  # 这里不同于kel.in中的 emax, 和scdft_tc, lambda_mu_kf中一样
            qe.write(" electron_maxstep = 100       \n")
            qe.write(" conv_thr = 1.0e-15           \n")
            qe.write(" spin_fluc =.true.            \n")  
            qe.write("/                             \n")      
        return inputfilename
    
    # not split mode
    def write_ph_no_split_in(self, work_directory:Path):
        inputfilename = "ph_no_split.in"
        ph_in = work_directory.joinpath(inputfilename)
        with open(ph_in, "w") as qe:
            qe.write("Electron-phonon coefficients for {}                \n".format(self.qe_inputpara.system_name))                                    
            qe.write(" &inputph                                          \n")      
            qe.write("  tr2_ph={},                                       \n".format(str(self.qe_inputpara.tr2_ph)))              
            qe.write("  prefix='{}',                                     \n".format(self.qe_inputpara.system_name))                
            # qe.write("  fildvscf='{}.dv',                                \n".format(self.qe_inputpara.system_name))       
            qe.write("  fildvscf='dvscf',                                \n")
            qe.write("  outdir='./tmp',                                  \n")               
            qe.write("  fildyn='{}.dyn',                                 \n".format(self.qe_inputpara.system_name))  
            qe.write("  alpha_mix(1)={},                                 \n".format(str(self.qe_inputpara.alpha_mix)))  # 可以修改的更小一些, 如果用vasp计算声子谱稳定, 可以修改为0.3
            qe.write("  ldisp=.true.,                                    \n") # true 代表通过 nq1 nq2 nq3方式计算声子，false 代表通过 坐标指定计算声子
            qe.write("  search_sym = {}                                  \n".format(self.qe_inputpara.search_sym))  # 如果是SCTK计算，最好用.false.         
            for i, species_name in enumerate(self.qe_inputpara.composition.keys()):
                element      = Element(species_name)
                species_mass = str(element.atomic_mass).strip("amu")
                qe.write("  amass({})={},                                \n".format(i+1, species_mass))             
            if self.qe_inputpara.SCTK_flag == True:
                qe.write("  lshift_q = .true.                            \n") # 移动q网格以避免Γ点的奇点。SCTK软件要求用四面体方法计算声子，这个是必须加的。
                if self.qe_inputpara.EPC_flag == True:
                    qe.write("  electron_phonon = 'scdft_input'          \n")
                    qe.write("  el_ph_sigma={},                          \n".format(str(self.qe_inputpara.el_ph_sigma)))                
                    qe.write("  el_ph_nsigma={},                         \n".format(str(self.qe_inputpara.el_ph_nsigma)))
                    qe.write("  trans={},                                \n".format(self.qe_inputpara.trans))
                    qe.write("  elph_nbnd_min = {}                       \n".format(self.qe_inputpara.elph_nbnd_min))
                    qe.write("  elph_nbnd_max = {}                       \n".format(self.qe_inputpara.elph_nbnd_max))
                    qe.write("  nk1={},nk2={},nk3={},                    \n".format(self.qe_inputpara.qpoints[0], self.qe_inputpara.qpoints[1], self.qe_inputpara.qpoints[2]))                 
                    qe.write("  nq1={},nq2={},nq3={},                    \n".format(self.qe_inputpara.qpoints[0], self.qe_inputpara.qpoints[1], self.qe_inputpara.qpoints[2]))                 
                else: # 这里代表只计算声子，不算电声耦合
                    qe.write("  nk1={},nk2={},nk3={},                    \n".format(self.qe_inputpara.kpoints_sparse[0], self.qe_inputpara.kpoints_sparse[1], self.qe_inputpara.kpoints_sparse[2]))                 
                    qe.write("  nq1={},nq2={},nq3={},                    \n".format(self.qe_inputpara.qpoints[0],        self.qe_inputpara.qpoints[1],        self.qe_inputpara.qpoints[2]))                 
            else: # 这里代表不使用SCTK的代码计算声子
                if self.qe_inputpara.EPC_flag == True:
                    qe.write("  electron_phonon='{}',                    \n".format(self.qe_inputpara.electron_phonon))                              
                    qe.write("  el_ph_sigma={},                          \n".format(str(self.qe_inputpara.el_ph_sigma)))                
                    qe.write("  el_ph_nsigma={},                         \n".format(str(self.qe_inputpara.el_ph_nsigma)))
                    qe.write("  trans={},                                \n".format(self.qe_inputpara.trans))
                    qe.write("  nq1={},nq2={},nq3={},                    \n".format(self.qe_inputpara.qpoints[0], self.qe_inputpara.qpoints[1], self.qe_inputpara.qpoints[2]))                 
                else:
                    qe.write("  nq1={},nq2={},nq3={},                    \n".format(self.qe_inputpara.qpoints[0], self.qe_inputpara.qpoints[1], self.qe_inputpara.qpoints[2]))                 
            qe.write("/                                                  \n")
        return inputfilename
    
    def write_dyn0(self, work_directory:Path):
        inputfilename = self.qe_inputpara.system_name+".dyn0"
        dyn0_path = work_directory.joinpath(inputfilename)
        with open(dyn0_path, "w") as qe:
            qe.write("{:<5} {:<5} {:<5}              \n".format(str(self.qe_inputpara.q1), str(self.qe_inputpara.q2), str(self.qe_inputpara.q3)))
            qe.write("{}                             \n".format(str(len(self.qe_inputpara.q_list))))
            for q in self.qe_inputpara.q_list:
                qe.write("{:<30}  {:<30}  {:<30}     \n".format(q[0], q[1], q[2]))
        return inputfilename

    # split mode1
    def write_split_ph_dyn0(self, work_directory:Path, q3):
        inputfilename = "split_ph.in"
        split_ph = work_directory.joinpath(inputfilename)
        with open(split_ph, "w") as qe:
            qe.write("Electron-phonon coefficients for {}                \n".format(self.qe_inputpara.system_name))                                    
            qe.write(" &inputph                                          \n")      
            qe.write("  tr2_ph={},                                       \n".format(str(self.qe_inputpara.tr2_ph)))              
            qe.write("  prefix='{}',                                     \n".format(self.qe_inputpara.system_name))                
            # qe.write("  fildvscf='{}.dv',                                \n".format(self.qe_inputpara.system_name))                     
            qe.write("  fildvscf='dvscf',                                \n")
            qe.write("  outdir='./tmp',                                  \n")               
            qe.write("  fildyn='{}.dyn',                                 \n".format(self.qe_inputpara.system_name))  
            qe.write("  alpha_mix(1)={},                                 \n".format(str(self.qe_inputpara.alpha_mix)))  # 可以修改的更小一些, 如果用vasp计算声子谱稳定, 可以修改为0.3
            qe.write("  ldisp=.false.,                                   \n") # true 代表通过 nq1 nq2 nq3方式计算声子，false 代表通过 坐标指定计算声子
            qe.write("  search_sym = {}                                  \n".format(self.qe_inputpara.search_sym))  # 如果是SCTK计算，最好用.false.         
            for i, species_name in enumerate(self.qe_inputpara.composition.keys()):
                element      = Element(species_name)
                species_mass = str(element.atomic_mass).strip("amu")
                qe.write("  amass({})={},                                \n".format(i+1, species_mass))
            if self.qe_inputpara.SCTK_flag == True:
                qe.write("  lshift_q = .true.                            \n") # 移动q网格以避免Γ点的奇点。SCTK软件要求用四面体方法计算声子，这个是必须加的。
                if self.qe_inputpara.EPC_flag == True:
                    qe.write("  electron_phonon = 'scdft_input'          \n")
                    qe.write("  el_ph_sigma={},                          \n".format(str(self.qe_inputpara.el_ph_sigma)))                
                    qe.write("  el_ph_nsigma={},                         \n".format(str(self.qe_inputpara.el_ph_nsigma)))
                    qe.write("  trans={},                                \n".format(self.qe_inputpara.trans))
                    qe.write("  elph_nbnd_min = {}                       \n".format(self.qe_inputpara.elph_nbnd_min))
                    qe.write("  elph_nbnd_max = {}                       \n".format(self.qe_inputpara.elph_nbnd_max))
                    qe.write("  nk1={},nk2={},nk3={},                    \n".format(self.qe_inputpara.kpoints_sparse[0], self.qe_inputpara.kpoints_sparse[1], self.qe_inputpara.kpoints_sparse[2]))                 
                    qe.write(" {:<30} {:<30} {:<30}                      \n".format(q3[0], q3[1], q3[2]))
                else: # 这里代表只计算声子，不算电声耦合
                    qe.write("  nk1={},nk2={},nk3={},                    \n".format(self.qe_inputpara.kpoints_sparse[0], self.qe_inputpara.kpoints_sparse[1], self.qe_inputpara.kpoints_sparse[2]))                 
                    qe.write(" {:<30} {:<30} {:<30}                      \n".format(q3[0], q3[1], q3[2]))
            else: # 这里代表不使用SCTK的代码计算声子
                if self.qe_inputpara.EPC_flag == True:
                    qe.write("  electron_phonon='{}',                    \n".format(self.qe_inputpara.electron_phonon))                              
                    qe.write("  el_ph_sigma={},                          \n".format(str(self.qe_inputpara.el_ph_sigma)))                
                    qe.write("  el_ph_nsigma={},                         \n".format(str(self.qe_inputpara.el_ph_nsigma)))
                    qe.write("  trans={},                                \n".format(self.qe_inputpara.trans))
                    qe.write(" {:<30} {:<30} {:<30}                      \n".format(q3[0], q3[1], q3[2]))
                else:
                    qe.write(" {:<30} {:<30} {:<30}                      \n".format(q3[0], q3[1], q3[2]))
            qe.write("/                                                  \n")

        return inputfilename

    # split mode2
    def write_split_phassignQ(self, work_directory:Path, start_q, last_q):
        inputfilename = "split_ph.in"
        split_ph_path = work_directory.joinpath(inputfilename)
        with open(split_ph_path, "w") as qe:        
            qe.write("Electron-phonon coefficients for {}                \n".format(self.qe_inputpara.system_name))                                    
            qe.write(" &inputph                                          \n")      
            qe.write("  tr2_ph={},                                       \n".format(str(self.qe_inputpara.tr2_ph)))              
            qe.write("  prefix='{}',                                     \n".format(self.qe_inputpara.system_name))                
            # qe.write("  fildvscf='{}.dv',                                \n".format(self.qe_inputpara.system_name))                     
            qe.write("  fildvscf='dvscf',                                \n")
            qe.write("  outdir='./tmp',                                  \n")               
            qe.write("  fildyn='{}.dyn',                                 \n".format(self.qe_inputpara.system_name))  
            qe.write("  alpha_mix(1)={},                                 \n".format(str(self.qe_inputpara.alpha_mix)))  # 可以修改的更小一些, 如果用vasp计算声子谱稳定, 可以修改为0.3
            qe.write("  ldisp=.true.,                                   \n") # true 代表通过 nq1 nq2 nq3方式计算声子，false 代表通过 坐标指定计算声子
            qe.write("  search_sym = {}                                  \n".format(self.qe_inputpara.search_sym))  # 如果是SCTK计算，最好用.false.         
            
            qe.write("  start_q={}                                       \n".format(start_q)) 
            qe.write("  last_q={}                                        \n".format(last_q))

            for i, species_name in enumerate(self.qe_inputpara.composition.keys()):
                element      = Element(species_name)
                species_mass = str(element.atomic_mass).strip("amu")
                qe.write("  amass({})={},                                \n".format(i+1, species_mass))
                
            if self.qe_inputpara.SCTK_flag == True:
                qe.write("  lshift_q = .true.                            \n") # 移动q网格以避免Γ点的奇点。SCTK软件要求用四面体方法计算声子，这个是必须加的。
                if self.qe_inputpara.EPC_flag == True:
                    qe.write("  electron_phonon = 'scdft_input'          \n")
                    qe.write("  el_ph_sigma={},                          \n".format(str(self.qe_inputpara.el_ph_sigma)))                
                    qe.write("  el_ph_nsigma={},                         \n".format(str(self.qe_inputpara.el_ph_nsigma)))
                    qe.write("  trans={},                                \n".format(self.qe_inputpara.trans))
                    qe.write("  elph_nbnd_min = {}                       \n".format(self.qe_inputpara.elph_nbnd_min))
                    qe.write("  elph_nbnd_max = {}                       \n".format(self.qe_inputpara.elph_nbnd_max))
                    qe.write("  nk1={},nk2={},nk3={},                    \n".format(self.qe_inputpara.qpoints[0], self.qe_inputpara.qpoints[1], self.qe_inputpara.qpoints[2]))                 
                    qe.write("  nq1={},nq2={},nq3={},                    \n".format(self.qe_inputpara.qpoints[0], self.qe_inputpara.qpoints[1], self.qe_inputpara.qpoints[2]))                 
                else: # 这里代表只计算声子，不算电声耦合
                    qe.write("  nk1={},nk2={},nk3={},                    \n".format(self.qe_inputpara.kpoints_sparse[0], self.qe_inputpara.kpoints_sparse[1], self.qe_inputpara.kpoints_sparse[2]))                 
                    qe.write("  nq1={},nq2={},nq3={},                    \n".format(self.qe_inputpara.qpoints[0],        self.qe_inputpara.qpoints[1],        self.qe_inputpara.qpoints[2]))                 
            else: # 这里代表不使用SCTK的代码计算声子
                if self.qe_inputpara.EPC_flag == True:
                    qe.write("  electron_phonon='{}',                    \n".format(self.qe_inputpara.electron_phonon))                              
                    qe.write("  el_ph_sigma={},                          \n".format(str(self.qe_inputpara.el_ph_sigma)))                
                    qe.write("  el_ph_nsigma={},                         \n".format(str(self.qe_inputpara.el_ph_nsigma)))
                    qe.write("  trans={},                                \n".format(self.qe_inputpara.trans))
                    qe.write("  nq1={},nq2={},nq3={},                    \n".format(self.qe_inputpara.qpoints[0], self.qe_inputpara.qpoints[1], self.qe_inputpara.qpoints[2]))                 
                else:
                    qe.write("  nq1={},nq2={},nq3={},                    \n".format(self.qe_inputpara.qpoints[0], self.qe_inputpara.qpoints[1], self.qe_inputpara.qpoints[2]))                 
            qe.write("/                                                  \n")

        return inputfilename

    def write_q2r_in(self, work_directory:Path):
        inputfilename = "q2r.in"
        q2r_in = work_directory.joinpath(inputfilename)
        with open(q2r_in, "w") as qe:
            qe.write("&input                      \n")             
            qe.write("  la2F = {},                \n".format(self.qe_inputpara.la2F))
            qe.write("  zasr = 'simple',          \n")                         
            qe.write("  fildyn = '{}.dyn'         \n".format(self.qe_inputpara.system_name))                               
            qe.write("  flfrc = '{}.fc',          \n".format(self.qe_inputpara.system_name))
            if self.qe_inputpara.EPC_flag == True:                    
                qe.write("  el_ph_nsigma={},      \n".format(str(self.qe_inputpara.el_ph_nsigma)))
            if self.qe_inputpara.SCTK_flag == True:
                qe.write("  lshift_q = .true.     \n") # 移动q网格以避免Γ点的奇点。SCTK软件要求用四面体方法计算声子，这个是必须加的。
            qe.write("/                           \n")        
        return inputfilename

    def write_matdyn_in(self, work_directory:Path):
        inputfilename = "matdyn.in"
        matdyn_in = work_directory.joinpath(inputfilename)
        special_qpoints_number  = len(self.qe_inputpara.path_name_coords)
        inserted_qpoints_number = self.qe_inputpara.qinserted
        with open(matdyn_in, "w") as qe:
            qe.write("&input                                             \n")               
            qe.write(" asr = 'simple',                                   \n")                        
            for i, species_name in enumerate(self.qe_inputpara.composition.keys()):
                element      = Element(species_name)
                species_mass = str(element.atomic_mass).strip("amu")
                qe.write(" amass({})={},                                 \n".format(i+1, species_mass))                                   
            qe.write(" flfrc = '{}.fc',                                 \n".format(self.qe_inputpara.system_name))                              
            qe.write(" flfrq='{}.freq',                                 \n".format(self.qe_inputpara.system_name))                              
            qe.write(" la2F = {},                                       \n".format(self.qe_inputpara.la2F))
            qe.write(" dos=.false.,                                     \n")                     
            qe.write(" q_in_band_form=.true.,                           \n")                                 
            qe.write(" q_in_cryst_coord=.true.,                         \n")
            if self.qe_inputpara.EPC_flag == True:                  
                qe.write(" el_ph_nsigma={},                                 \n".format(str(self.qe_inputpara.el_ph_nsigma)))
            if self.qe_inputpara.SCTK_flag == True:
                qe.write("  lshift_q = .true.                            \n") # 移动q网格以避免Γ点的奇点。SCTK软件要求用四面体方法计算声子，这个是必须加的。
            qe.write("/                                                  \n")          
            qe.write("{}                                                 \n".format(special_qpoints_number))            
            for name, coord in self.qe_inputpara.path_name_coords:
                qe.write(" {:<15} {:<15} {:<15} {:<5}                   \n".format(str(coord[0]), str(coord[1]), str(coord[2]), str(inserted_qpoints_number)))
        return inputfilename

    def write_phonobanddata_in(self, work_directory:Path):
        inputfilename = "phonobanddata.in"
        with open(inputfilename, "w") as qe:
            qe.write("{}.freq                    \n".format(self.qe_inputpara.system_name))            
            qe.write("0 4000                     \n".format())    
            qe.write("{}.phonon.bands.dat        \n".format(self.qe_inputpara.system_name))                        
            qe.write("{}.phonon.bands.ps         \n".format(self.qe_inputpara.system_name))                       
            qe.write("0                          \n".format())
            qe.write("50 0                       \n".format())  
        return inputfilename
    
    def write_eletdos_in(self, work_directory:Path):
        inputfilename = "eletdos.in"
        eledos_in = work_directory.joinpath(inputfilename) 
        with open(eledos_in, "w") as qe:
            qe.write("&dos                       \n")
            qe.write("   prefix = '{}',          \n".format(self.qe_inputpara.system_name))                                                          
            qe.write("   outdir = './tmp',       \n")                                                                                    
            qe.write("   fildos = '{}.tdos',     \n".format(self.qe_inputpara.system_name))
            qe.write("   DeltaE = {},            \n".format(self.qe_inputpara.DeltaE))                                    
            qe.write("   emin = {},              \n".format(self.qe_inputpara.emin))
            qe.write("   emax = {},              \n".format(self.qe_inputpara.emax))
            qe.write("/                          \n")                                                                
        return inputfilename

    def write_elepdos_in(self, work_directory:Path):
        inputfilename = "elepdos.in"
        eledos_in = work_directory.joinpath(inputfilename) 
        with open(eledos_in, "w") as qe:
            qe.write("&projwfc                   \n")
            qe.write("   prefix = '{}',          \n".format(self.qe_inputpara.system_name))                                                          
            qe.write("   outdir = './tmp',       \n")                                                                                    
            qe.write("   filpdos= '{}.pdos',     \n".format(self.qe_inputpara.system_name))
            qe.write("   ngauss = {}             \n".format(self.qe_inputpara.ngauss))
            qe.write("   filproj= '{}.proj',     \n".format(self.qe_inputpara.system_name))
            qe.write("   DeltaE = {},            \n".format(self.qe_inputpara.DeltaE))   
            qe.write("   degauss = {},           \n".format(self.qe_inputpara.pdosdegauss))                                 
            qe.write("   emin = {},              \n".format(self.qe_inputpara.emin))
            qe.write("   emax = {},              \n".format(self.qe_inputpara.emax))
            qe.write("/                          \n")                                                                
        return inputfilename

    def write_eleband_in(self, work_directory:Path):
        inputfilename = "eleband.in"
        eleband_in = work_directory.joinpath(inputfilename)
        inserted_kpoints_number = self.qe_inputpara.kinserted
        with open(eleband_in, "w") as qe:
            qe.write("&CONTROL\n")
            qe.write(" calculation='bands',            \n")
            qe.write(" restart_mode='from_scratch',    \n")
            qe.write(" prefix='{}',                    \n".format(self.qe_inputpara.system_name))
            qe.write(" pseudo_dir='{}',                \n".format(str(self.qe_inputpara.workpath_pppath.absolute())))
            qe.write(" verbosity = 'high',             \n")  
            qe.write(" outdir='./tmp',                 \n")
            qe.write("/\n")

            qe.write("&SYSTEM\n")
            qe.write(" ibrav=0,                        \n")  # 设置ibrav=0，这时需要在输入文件中写入CELL_PARAMETERS，即CELL的基矢量. alat bohr angstrom alat 由 celldm(1)或A定义的晶格常数单位
            qe.write(" nat={},                         \n".format(self.qe_inputpara.all_atoms_quantity))
            qe.write(" ntyp={},                        \n".format(self.qe_inputpara.species_quantity))
            qe.write(" occupations = '{}',             \n".format(self.qe_inputpara.occupations))
            if self.qe_inputpara.occupations == 'smearing':
                qe.write(" occupations = '{}',         \n".format(self.qe_inputpara.occupations)) 
                qe.write(" smearing = '{}',            \n".format(self.qe_inputpara.smearing))
                qe.write(" degauss = {},               \n".format(self.qe_inputpara.degauss))
            elif self.qe_inputpara.occupations == 'tetrahedra' or self.qe_inputpara.occupations == 'tetrahedra_opt':
                qe.write(" occupations = '{}',         \n".format(self.qe_inputpara.occupations)) 
            qe.write(" ecutwfc = {},                   \n".format(self.qe_inputpara.ecutwfc))
            qe.write(" ecutrho = {},                   \n".format(self.qe_inputpara.ecutrho))
            qe.write(" lspinorb = .{}.,                \n".format(self.qe_inputpara.lspinorb))
            qe.write(" noncolin = .{}.,                \n".format(self.qe_inputpara.noncolin))
            if self.qe_inputpara.celldm1 is not None:
                qe.write(" celldm(1) = {:<.16f},       \n".format(self.qe_inputpara.celldm1))
            if self.qe_inputpara.nbnd is not None:
                qe.write(" nbnd = {},                  \n".format(self.qe_inputpara.nbnd))
            qe.write("/\n")

            qe.write("&ELECTRONS\n")
            qe.write(" conv_thr = {},                  \n".format(self.qe_inputpara.conv_thr))
            qe.write(" mixing_beta = {},               \n".format(self.qe_inputpara.mixing_beta))
            qe.write(" electron_maxstep = {},          \n".format(self.qe_inputpara.electron_maxstep))
            qe.write("/\n")

            qe.write("ATOMIC_SPECIES                   \n")
            for species_name in self.qe_inputpara.composition.keys():
                species_pseudo = get_pps_for_a_element(species_name, self.qe_inputpara.final_choosed_pp)
                if len(species_pseudo)==1:
                    print(f"write USPP for species in scf.in: {species_pseudo[0]}") 
                    element      = Element(species_name)
                    species_mass = str(element.atomic_mass).strip("amu")
                    qe.write(" {:<5}  {:<10}  {:<50} \n".format(species_name, species_mass, species_pseudo[0]))
                else:
                    print("\nNote: --------------------")
                    print("    When write pseudo-potentional information, something wrong. Maybe methdo `get_pps_for_a_element` is problematic !")
                    sys.exit(1)
            if "V3_Hessian" in self.qe_inputpara.input_file_path.name:
                qe.write("CELL_PARAMETERS {alat}          \n")  # 如果选择alat单位，原子坐标选择分数坐标，即，ATOMIC_POSITIONS (crystal), 且不设置celldm(1). 这时alat和celldm(1)设置成v1的长度
            else:
                qe.write("CELL_PARAMETERS {angstrom}      \n")  # 如果选择angstrom单位，原子坐标选择分数坐标，即，ATOMIC_POSITIONS (crystal), 且不设置celldm(1). 这时alat和celldm(1)设置成v1的长度
            for cell_p in self.qe_inputpara.cell_parameters:
                qe.write("{}\n".format(cell_p))
            if "V3_Hessian" in self.qe_inputpara.input_file_path.name:
                qe.write("ATOMIC_POSITIONS {alat}         \n")
            else:
                qe.write("ATOMIC_POSITIONS {crystal}      \n")  # 如果选择angstrom单位，原子坐标选择分数坐标，即，ATOMIC_POSITIONS (crystal), 且不设置celldm(1). 这时alat和celldm(1)设置成v1的长度
            for site in self.qe_inputpara.fractional_sites:
                qe.write("{}\n".format(site))
            qe.write("K_POINTS {crystal_b}             \n")
            qe.write(" {}                              \n".format(len(self.qe_inputpara.path_name_coords)))
            for name, coord in self.qe_inputpara.path_name_coords:
                qe.write(" {:<15} {:<15} {:<15}  {}  ! {:<5}                   \n".format(str(coord[0]), str(coord[1]), str(coord[2]), str(inserted_kpoints_number), name))
        return inputfilename

    def write_elebanddata_in(self, work_directory:Path):
        inputfilename = "elebanddata.in"
        eletronband_in = work_directory.joinpath(inputfilename)
        with open(eletronband_in, "w") as f:
            f.write("&BANDS\n")
            f.write(" prefix='{}',\n".format(self.qe_inputpara.system_name))
            f.write(" outdir='./tmp',\n")
            f.write(" filband='elebanddata.dat',\n")
            f.write(" lp=.true.\n")
            f.write("/\n")
        return inputfilename
    
    def write_elebandprojdata_in(self, work_directory:Path):
        inputfilename = "elebandprojdata.in"
        eletronband_in = work_directory.joinpath(inputfilename)
        with open(eletronband_in, "w") as f:
            f.write("&projwfc\n")
            f.write(" prefix='{}',\n".format(self.qe_inputpara.system_name))
            f.write(" outdir='./tmp',\n")
            f.write(" lsym=.{}.,\n".format(self.qe_inputpara.lsym))
            f.write(" filproj = 'elebandprojdata'\n")
            f.write("/\n")
        return inputfilename

    def write_phonodos_in(self, work_directory:Path):
        inputfilename = "phonodos.in"
        phonodos_in = work_directory.joinpath(inputfilename) 
        with open(phonodos_in, "w") as qe:
            qe.write("&input                      \n")
            qe.write("   asr = 'simple',          \n")                                 
            for i, species_name in enumerate(self.qe_inputpara.composition.keys()):
                element      = Element(species_name)
                species_mass = str(element.atomic_mass).strip("amu")
                qe.write("   amass({})={},                               \n".format(i+1, species_mass))                                    
            qe.write("   flfrc = '{}.fc',         \n".format(self.qe_inputpara.system_name))
            qe.write("   flfrq = '{}.freq',       \n".format(self.qe_inputpara.system_name))                                         
            qe.write("   la2F = {},               \n".format(self.qe_inputpara.la2F))
            qe.write("   dos = .true.,            \n")  # 计算声子态密度,dos必须设置为.true.                           
            qe.write("   fldos = '{}.dos',        \n".format(self.qe_inputpara.system_name+"_phono"))                                       
            qe.write("   nk1={}, nk2={}, nk3={},  \n".format(self.qe_inputpara.qpoints[0], self.qe_inputpara.qpoints[1], self.qe_inputpara.qpoints[2]))  # 计算态密度时要用更密的q点网格，这需设置nk1, nk2, nk3
            qe.write("   ndos={},                 \n".format(self.qe_inputpara.ndos))  # 态密度的能量刻度上的点的数目                       
            qe.write("   el_ph_nsigma={},         \n".format(self.qe_inputpara.el_ph_nsigma))
            qe.write("/                           \n")                                                                
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
            top_freq        = self.qe_inputpara.top_freq
            broaden         = self.qe_inputpara.broaden
            smearing_method = self.qe_inputpara.smearing_method
            print("\nNote: --------------------")
            print("    Again, check to see if the four values are the same")
            print("    qirreduced_coords:{}   qweights:{}   qirreduced number:{}   elphInpLambda number:{}".format(len(self.qe_inputpara.qirreduced_coords),  len(self.qe_inputpara.qweights), int(self.qe_inputpara.qirreduced), len(elphInpLambda)))
            # time.sleep(3)
            if len(self.qe_inputpara.qirreduced_coords) == len(self.qe_inputpara.qweights) == int(self.qe_inputpara.qirreduced) == len(elphInpLambda):
                q_number = self.qe_inputpara.qirreduced
                q_coords = self.qe_inputpara.qirreduced_coords
                q_weight = self.qe_inputpara.qweights
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

        temperature_steps = self.qe_inputpara.temperature_steps
        inputfilename = "INPUT"
        eliashberg_in = work_directory.joinpath(inputfilename)
        with open(eliashberg_in, "w") as qe:
            qe.write("{:<10} {:<10}".format(screen_constant, temperature_steps))
        return inputfilename

    def write_alpha2f_out(self, work_directory:Path, gaussid):
        alpha2f_out = work_directory.joinpath("ALPHA2F.OUT").absolute()
        alpha2F_dat = work_directory.joinpath("alpha2F.dat").absolute()
        
        # 方法1 通过计算声子态密度, 获得a2F.dos*文件, 然后选择一个合适的文件用来计算超导
        if self.qe_inputpara.a2fdos:
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

    def write_epw_energyband_in(self, work_directory:Path):
        inputfilename = "epw_energyband.in"
        epw_energyband_in = work_directory.joinpath(inputfilename)
        with open(epw_energyband_in, "w") as epw:
            epw.write("&inputepw\n")
            epw.write(" prefix      ='{}',\n".format(self.qe_inputpara.system_name))
            for i, species_name in enumerate(self.qe_inputpara.composition.keys()):
                element      = Element(species_name)
                species_mass = str(element.atomic_mass).strip("amu")
                epw.write(" amass({})    ={},\n".format(i+1, species_mass))          
            epw.write(" outdir      ='./tmp'\n")
            epw.write(" dvscf_dir   ='{}'\n".format(self.qe_inputpara.dvscf_dir))

            epw.write(" etf_mem     ={}\n".format(self.qe_inputpara.etf_mem))

            epw.write(" use_ws      = .false.\n")
            epw.write(" wannierize  = .true.\n")
            epw.write(" nbndsub     = {}\n".format(self.qe_inputpara.nbndsub))
            epw.write(" bands_skipped = 'exclude_bands = {}'\n".format(self.qe_inputpara.bands_skipped))
            epw.write(" num_iter    = {}\n".format(self.qe_inputpara.num_iter))
            epw.write(" dis_froz_min = {}\n".format(self.qe_inputpara.dis_froz_min))
            epw.write(" dis_froz_max = {}\n".format(self.qe_inputpara.dis_froz_max))
            for idx, pj in enumerate(self.qe_inputpara.proj):
                epw.write(" proj({})     = {}\n".format(idx+1, pj))
            epw.write(" wannier_plot= .true.\n")
            epw.write(" wdata(1)    = 'bands_plot = .true.'\n")
            epw.write(" wdata(2)    = 'begin kpoint_path'\n")
            for idx, path_name_coord in enumerate(self.qe_inputpara.path_name_coords_for_EPW):
                epw.write(f" wdata({idx+3})    = '{path_name_coord}'\n")
            epw.write(f" wdata({idx+4})    = 'end kpoint_path'\n")
            epw.write("/                           \n")                                                                
        
        return inputfilename
    
    def write_epw_phono_in(self, work_directory:Path):
        inputfilename = "epw_phono.in"
        epw_phono_in = work_directory.joinpath(inputfilename)
        with open(epw_phono_in, "w") as epw:
            epw.write("&inputepw\n")
            epw.write(" prefix      ='{}',\n".format(self.qe_inputpara.system_name))
            for i, species_name in enumerate(self.qe_inputpara.composition.keys()):
                element      = Element(species_name)
                species_mass = str(element.atomic_mass).strip("amu")
                epw.write(" amass({})    ={},\n".format(i+1, species_mass))          
            epw.write(" outdir      ='./tmp'\n")
            epw.write(" dvscf_dir   ='{}'".format(self.qe_inputpara.dvscf_dir))

            epw.write(" etf_mem     ={}\n".format(self.qe_inputpara.etf_mem))
            epw.write(" elph        = .true.\n")
            epw.write(" epbwrite    = .true.\n")
            epw.write(" epbread     = .false.\n")
            epw.write(" epwwrite    = .true.\n")
            epw.write(" epwread     = .false.\n")

            epw.write(" use_ws      = .false.\n")
            epw.write(" wannierize  = .true.\n")
            epw.write(" nbndsub     = {}\n".format(self.qe_inputpara.nbndsub))
            epw.write(" bands_skipped = 'exclude_bands = {}'\n".format(self.qe_inputpara.bands_skipped))
            epw.write(" num_iter    = {}\n".format(self.qe_inputpara.num_iter))
            epw.write(" dis_froz_min = {}\n".format(self.qe_inputpara.dis_froz_min))
            epw.write(" dis_froz_max = {}\n".format(self.qe_inputpara.dis_froz_max))
            for idx, pj in enumerate(self.qe_inputpara.proj):
                epw.write(" proj({})     = {}\n".format(idx, pj))
            epw.write(" wannier_plot= .true.\n")
            epw.write(" wdata(1)    = 'bands_plot = .true.'\n")
            epw.write(" wdata(2)    = 'begin kpoint_path'\n")
            for idx, path_name_coord in enumerate(self.qe_inputpara.path_name_coords_for_EPW):
                epw.write(f" wdata({idx+3})    = '{path_name_coord}'\n")
            epw.write(f" wdata({idx+4})    = 'end kpoint_path'\n")
        return inputfilename

    def write_epw_elph_in():
        pass

    def write_epw_aniso_sc_in():
        pass
