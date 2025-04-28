import os
import re
import sys
import time
import logging
from pathlib import Path

from epw.epw_inputpara import epw_inputpara

from pymatgen.core.periodic_table import Element

logger = logging.getLogger("epw_writeinput")

class epw_writeinput:
    
    def __init__(
        self, 
        epw_inputpara: epw_inputpara
        ):

        self.epw_inputpara = epw_inputpara

    # write inputfile
    def writeinput(self, mode=None):
        if mode == None:
            mode = self.epw_inputpara.mode

        if mode == "epw_eband":
            inputfilename = self.write_epw_eband_in(self.epw_inputpara.work_path)
            return inputfilename
        if mode == "epw_phono":
            inputfilename1 = self.write_epw_phonodata_in(self.epw_inputpara.work_path)
            inputfilename2 = self.write_epw_phonodata_plot_in(self.epw_inputpara.work_path)
            return inputfilename1, inputfilename2
        if mode == "epw_elph":
            inputfilename = self.write_epw_elph_in(self.epw_inputpara.work_path)
            return inputfilename
        if mode == "epw_sc":
            inputfilename1 = self.write_epw_iso_sc_in(self.epw_inputpara.work_path)
            inputfilename2 = self.write_epw_aniso_sc_in(self.epw_inputpara.work_path)
            return inputfilename1, inputfilename2

    def write_epw_eband_in(self, work_directory:Path):
        inputfilename = "epw_eband.in"
        epw_energyband_in = work_directory.joinpath(inputfilename)
        with open(epw_energyband_in, "w") as epw:
            epw.write("&inputepw\n")
            epw.write(" prefix      ='{}',\n".format(self.epw_inputpara.system_name))
            for i, species_name in enumerate(self.epw_inputpara.composition.keys()):
                element      = Element(species_name)
                species_mass = str(element.atomic_mass).strip("amu")
                epw.write(" amass({})    ={},\n".format(i+1, species_mass))          
            epw.write(" outdir      ='./tmp'\n")
            epw.write(" dvscf_dir   ='{}'\n".format(self.epw_inputpara.dvscf_dir))
            epw.write("\n")
            epw.write(" etf_mem     ={}\n".format(self.epw_inputpara.etf_mem))
            epw.write("\n")
            epw.write(" elph        = .false.\n")
            epw.write("\n")
            epw.write(" use_ws      = .true.\n")
            epw.write(" wannierize  = .true.\n")
            epw.write(" nbndsub     = {}\n".format(self.epw_inputpara.nbndsub))
            if self.epw_inputpara.bands_skipped is not None:
                epw.write(" bands_skipped = 'exclude_bands = {}'\n".format(self.epw_inputpara.bands_skipped))
            epw.write(" num_iter    = {}\n".format(self.epw_inputpara.num_iter))
            epw.write("\n")
            epw.write(" dis_win_max = {}\n".format(self.epw_inputpara.dis_win_max))
            epw.write(" dis_froz_min = {}\n".format(self.epw_inputpara.dis_froz_min))
            epw.write(" dis_froz_max = {}\n".format(self.epw_inputpara.dis_froz_max))
            epw.write("\n")
            for idx, pj in enumerate(self.epw_inputpara.proj):
                epw.write(" proj({})     = {}\n".format(idx+1, pj))
            epw.write(" wannier_plot= .true.\n")
            epw.write(" wdata(1)    = 'bands_plot = .true.'\n")
            epw.write(" wdata(2)    = 'begin kpoint_path'\n")
            for idx, path_name_coord in enumerate(self.epw_inputpara.path_name_coords_for_EPW):
                epw.write(f" wdata({idx+3})    = '{path_name_coord}'\n")
            epw.write(f" wdata({idx+4})    = 'end kpoint_path'\n")
            epw.write("\n")
            epw.write(" wdata({})    = 'dis_num_iter = {}'\n".format(idx+5, self.epw_inputpara.dis_num_iter))
            epw.write("\n")
            epw.write(" nk1 = {}\n nk2 = {}\n nk3 = {}\n".format(self.epw_inputpara.nk[0], self.epw_inputpara.nk[1], self.epw_inputpara.nk[2]))
            epw.write("/                           \n")                                                                
        
        return inputfilename
    
    def write_epw_phono_in(self, work_directory:Path):
        inputfilename = "epw_phono.in"
        epw_phono_in = work_directory.joinpath(inputfilename)
        with open(epw_phono_in, "w") as epw:
            
            epw.write("&inputepw\n")
            epw.write(" prefix      ='{}',\n".format(self.epw_inputpara.system_name))
            for i, species_name in enumerate(self.epw_inputpara.composition.keys()):
                element      = Element(species_name)
                species_mass = str(element.atomic_mass).strip("amu")
                epw.write(" amass({})    ={},\n".format(i+1, species_mass))          
            epw.write(" outdir      ='./tmp'\n")
            epw.write(" dvscf_dir   ='{}'\n".format(self.epw_inputpara.dvscf_dir))
            epw.write("\n")
            epw.write(" etf_mem     ={}\n".format(self.epw_inputpara.etf_mem))
            epw.write("\n")
            epw.write(" elph        = .true.\n")
            epw.write(" epbwrite    = .true.\n")
            epw.write(" epbread     = .false.\n")
            epw.write(" epwwrite    = .true.\n")
            epw.write(" epwread     = .false.\n")
            epw.write("\n")
            epw.write(" use_ws      = .true.\n")
            epw.write(" wannierize  = {}\n".format(self.epw_inputpara.wannierize))
            epw.write(" nbndsub     = {}\n".format(self.epw_inputpara.nbndsub))
            if self.epw_inputpara.bands_skipped is not None:
                epw.write(" bands_skipped = 'exclude_bands = {}'\n".format(self.epw_inputpara.bands_skipped))
            epw.write(" num_iter    = {}\n".format(self.epw_inputpara.num_iter))
            epw.write("\n")
            epw.write(" dis_win_max = {}\n".format(self.epw_inputpara.dis_win_max))
            epw.write(" dis_froz_min = {}\n".format(self.epw_inputpara.dis_froz_min))
            epw.write(" dis_froz_max = {}\n".format(self.epw_inputpara.dis_froz_max))
            epw.write("\n")
            for idx, pj in enumerate(self.epw_inputpara.proj):
                epw.write(" proj({})     = {}\n".format(idx+1, pj))
                
            epw.write(" wdata({})    = 'dis_num_iter = {}'\n".format(idx+1, self.epw_inputpara.dis_num_iter))
            epw.write("\n")
            epw.write(" nk1 = {}\n nk2 = {}\n nk3 = {}\n".format(self.epw_inputpara.nk[0], self.epw_inputpara.nk[1], self.epw_inputpara.nk[2]))
            epw.write("\n")
            epw.write(" nq1 = {}\n nq2 = {}\n nq3 = {}\n".format(self.epw_inputpara.nq[0], self.epw_inputpara.nq[1], self.epw_inputpara.nq[2]))
            epw.write("\n")
            epw.write(" nkf1 = {}\n nkf2 = {}\n nkf3 = {}\n".format(1,1,1))
            epw.write("\n")
            epw.write(" nqf1 = {}\n nqf2 = {}\n nqf3 = {}\n".format(1,1,1))

            epw.write("/                           \n")                                                                

        return inputfilename

    def write_epw_phonodata_in(self, work_directory:Path):
        inputfilename = "epw_phonodata.in"
        epw_phonodata_in = work_directory.joinpath(inputfilename)
        with open(epw_phonodata_in, "w") as epw:
            
            epw.write("&inputepw\n")
            epw.write(" prefix      ='{}',\n".format(self.epw_inputpara.system_name))
            for i, species_name in enumerate(self.epw_inputpara.composition.keys()):
                element      = Element(species_name)
                species_mass = str(element.atomic_mass).strip("amu")
                epw.write(" amass({})    ={},\n".format(i+1, species_mass))          
            epw.write(" outdir      ='./tmp'\n")
            epw.write(" dvscf_dir   ='{}'\n".format(self.epw_inputpara.dvscf_dir))
            epw.write("\n")
            epw.write(" etf_mem     ={}\n".format(self.epw_inputpara.etf_mem))
            epw.write("\n")
            epw.write(" elph        = .true.\n")
            epw.write(" epbwrite    = .false.\n")
            epw.write(" epbread     = .false.\n")
            epw.write(" epwwrite    = .false.\n")
            epw.write(" epwread     = .true.\n")
            epw.write("\n")
            epw.write(" lifc = {}\n".format(self.epw_inputpara.lifc))
            epw.write("\n")
            epw.write(" use_ws      = .true.\n")
            epw.write(" wannierize  = {}\n".format(self.epw_inputpara.wannierize))
            epw.write(" nbndsub     = {}\n".format(self.epw_inputpara.nbndsub))
            if self.epw_inputpara.bands_skipped is not None:
                epw.write(" bands_skipped = 'exclude_bands = {}'\n".format(self.epw_inputpara.bands_skipped))
            epw.write(" num_iter    = {}\n".format(self.epw_inputpara.num_iter))
            epw.write("\n")
            epw.write(" dis_win_max = {} \n".format(self.epw_inputpara.dis_win_max))
            epw.write(" dis_froz_min = {}\n".format(self.epw_inputpara.dis_froz_min))
            epw.write(" dis_froz_max = {}\n".format(self.epw_inputpara.dis_froz_max))
            epw.write("\n")
            for idx, pj in enumerate(self.epw_inputpara.proj):
                epw.write(" proj({})     = {}\n".format(idx+1, pj))
                
            epw.write(" wdata({})    = 'dis_num_iter = {}'\n".format(idx+1, self.epw_inputpara.dis_num_iter))
            epw.write(" band_plot = .true. \n")
            epw.write(" filkf = '{}_band.kpt'\n".format(self.epw_inputpara.system_name))
            epw.write(" filqf = '{}_band.kpt'\n".format(self.epw_inputpara.system_name))
            epw.write(" asr_typ = '{}'     \n".format(self.epw_inputpara.asr_typ))

            epw.write("\n")
            epw.write(" nk1 = {}\n nk2 = {}\n nk3 = {}\n".format(self.epw_inputpara.nk[0], self.epw_inputpara.nk[1], self.epw_inputpara.nk[2]))
            epw.write("\n")
            epw.write(" nq1 = {}\n nq2 = {}\n nq3 = {}\n".format(self.epw_inputpara.nq[0], self.epw_inputpara.nq[1], self.epw_inputpara.nq[2]))
            epw.write("/                           \n")                                                                

        return inputfilename

    def write_epw_phonodata_plot_in(self, work_directory:Path):
        inputfilename = "epw_phonodata_plot.in"
        epw_phonodata_plot_in = work_directory.joinpath(inputfilename)
        with open(epw_phonodata_plot_in, "w") as epw:
            epw.write("phband.freq\n")
            epw.write("0 500\n")
            epw.write("freq.dat\n")
            epw.write("\n")
            epw.write("# Second line is min and max value in your spectrum\n")
            epw.write("# Meanwhile, another file named band.eig is eletron spectrum.\n")
            epw.write("# band.eig also can be plotted by plotband.x, detailed in Tue.6.Lafuente.pdf.\n") 
        return inputfilename

    def write_epw_elph_in(self, work_directory:Path):
        inputfilename = "epw_elph.in"
        epw_phonodata_in = work_directory.joinpath(inputfilename)
        with open(epw_phonodata_in, "w") as epw:
            
            epw.write("&inputepw\n")
            epw.write(" prefix      ='{}',\n".format(self.epw_inputpara.system_name))
            for i, species_name in enumerate(self.epw_inputpara.composition.keys()):
                element      = Element(species_name)
                species_mass = str(element.atomic_mass).strip("amu")
                epw.write(" amass({})    ={},\n".format(i+1, species_mass))          
            epw.write(" outdir      ='./tmp'\n")
            epw.write(" dvscf_dir   ='{}'\n".format(self.epw_inputpara.dvscf_dir))
            epw.write("\n")
            epw.write(" etf_mem     ={}\n".format(self.epw_inputpara.etf_mem))
            epw.write("\n")
            epw.write(" ep_coupling = .true.\n")
            epw.write(" elph        = .true.\n")
            epw.write(" epbwrite    = .false.  ! electron-phonon matrix elements in the coarse Bloch representation and relevant data (dyn matrices)\n")
            epw.write(" epbread     = .false.  ! electron-phonon matrix elements in the coarse Bloch representation and relevant data (dyn matrices)\n")
            epw.write(" epwwrite    = .false.  ! electron-phonon matrix elements in the coarse Wannier representation and relevant data (dyn matrices)\n")
            epw.write(" epwread     = .true.   ! electron-phonon matrix elements in the coarse Wannier representation and relevant data (dyn matrices). ! It is used for a restart calculation and requires kmaps = .true.  prefix.epmatwp, crystal.fmt, vmedata.fmt, wigner.fmt, and epwdata.fmt are all needed for restarting a calculation.\n")
            epw.write("\n")
            epw.write("\n")
            epw.write(" lifc = {}\n".format(self.epw_inputpara.lifc))
            epw.write("\n")
            epw.write(" lifc = {}\n".format(self.epw_inputpara.lifc))
            epw.write(" wannierize  = {}\n".format(self.epw_inputpara.wannierize))
            epw.write(" nbndsub     = {}\n".format(self.epw_inputpara.nbndsub))
            if self.epw_inputpara.bands_skipped is not None:
                epw.write(" bands_skipped = 'exclude_bands = {}'\n".format(self.epw_inputpara.bands_skipped))
            epw.write(" num_iter    = {}\n".format(self.epw_inputpara.num_iter))
            epw.write("\n")
            epw.write(" dis_win_max = {} \n".format(self.epw_inputpara.dis_win_max))
            epw.write(" dis_froz_min = {}\n".format(self.epw_inputpara.dis_froz_min))
            epw.write(" dis_froz_max = {}\n".format(self.epw_inputpara.dis_froz_max))
            epw.write("\n")
            for idx, pj in enumerate(self.epw_inputpara.proj):
                epw.write(" proj({})     = {}\n".format(idx+1, pj))
                
            epw.write(" wdata({})    = 'dis_num_iter = {}'\n".format(idx+1, self.epw_inputpara.dis_num_iter))
            epw.write("\n")
            
            epw.write(" iverbosity  = 2     ! 2 = verbose output for the superconducting part only.\n")
            epw.write(" fsthick  = {}       ! Fermi window thickness [eV]\n".format(self.epw_inputpara.fsthick))
            epw.write(" degaussw = {}       ! smearing in energy-conserving delta functions in [eV]\n".format(self.epw_inputpara.degaussw))
            epw.write(" degaussq = {}       ! smearing for sum over q in the e-ph coupling in [meV]\n".format(self.epw_inputpara.degaussq))
            epw.write(" ephwrite  = .true.  ! write ephmatXX, egnv, freq, and ikmap files in prefix.ephmat directory\n")

            epw.write("\n")
            epw.write(" nk1 = {}\n nk2 = {}\n nk3 = {}\n".format(self.epw_inputpara.nk[0], self.epw_inputpara.nk[1], self.epw_inputpara.nk[2]))
            epw.write("\n")
            epw.write(" nq1 = {}\n nq2 = {}\n nq3 = {}\n".format(self.epw_inputpara.nq[0], self.epw_inputpara.nq[1], self.epw_inputpara.nq[2]))
            epw.write("\n")
            epw.write(" nkf1 = {}\n nkf2 = {}\n nkf3 = {}\n".format(self.epw_inputpara.nkf[0], self.epw_inputpara.nkf[1], self.epw_inputpara.nkf[2]))
            epw.write("\n")
            epw.write(" nqf1 = {}\n nqf2 = {}\n nqf3 = {}\n".format(self.epw_inputpara.nqf[0], self.epw_inputpara.nqf[1], self.epw_inputpara.nqf[2]))
            epw.write("/                           \n")                                                                

        return inputfilename

    def write_epw_iso_sc_in(self, work_directory:Path):
        inputfilename = "epw_iso_sc.in"
        epw_phonodata_in = work_directory.joinpath(inputfilename)
        with open(epw_phonodata_in, "w") as epw:
            
            epw.write("&inputepw\n")
            epw.write(" prefix      ='{}',\n".format(self.epw_inputpara.system_name))
            for i, species_name in enumerate(self.epw_inputpara.composition.keys()):
                element      = Element(species_name)
                species_mass = str(element.atomic_mass).strip("amu")
                epw.write(" amass({})    ={},\n".format(i+1, species_mass))          
            epw.write(" outdir      ='./tmp'\n")
            epw.write(" dvscf_dir   ='{}'\n".format(self.epw_inputpara.dvscf_dir))
            epw.write("\n")
            epw.write(" etf_mem     ={}\n".format(self.epw_inputpara.etf_mem))
            epw.write("\n")
            epw.write(" ep_coupling = .false.\n")
            epw.write(" elph        = .false.\n")
            epw.write(" epbwrite    = .false.  ! electron-phonon matrix elements in the coarse Bloch representation and relevant data (dyn matrices)\n")
            epw.write(" epbread     = .false.  ! electron-phonon matrix elements in the coarse Bloch representation and relevant data (dyn matrices)\n")
            epw.write(" epwwrite    = .false.  ! electron-phonon matrix elements in the coarse Wannier representation and relevant data (dyn matrices)\n")
            epw.write(" epwread     = .true.   ! electron-phonon matrix elements in the coarse Wannier representation and relevant data (dyn matrices). ! It is used for a restart calculation and requires kmaps = .true.  prefix.epmatwp, crystal.fmt, vmedata.fmt, wigner.fmt, and epwdata.fmt are all needed for restarting a calculation.\n")
            epw.write("\n")
            epw.write("\n")
            epw.write(" lifc = {}\n".format(self.epw_inputpara.lifc))
            epw.write("\n")
            epw.write(" use_ws      = .true.\n")
            epw.write(" wannierize  = {}\n".format(self.epw_inputpara.wannierize))
            epw.write(" nbndsub     = {}\n".format(self.epw_inputpara.nbndsub))
            if self.epw_inputpara.bands_skipped is not None:
                epw.write(" bands_skipped = 'exclude_bands = {}'\n".format(self.epw_inputpara.bands_skipped))
            epw.write(" num_iter    = {}\n".format(self.epw_inputpara.num_iter))
            epw.write("\n")
            epw.write(" dis_win_max = {} \n".format(self.epw_inputpara.dis_win_max))
            epw.write(" dis_froz_min = {}\n".format(self.epw_inputpara.dis_froz_min))
            epw.write(" dis_froz_max = {}\n".format(self.epw_inputpara.dis_froz_max))
            epw.write("\n")
            for idx, pj in enumerate(self.epw_inputpara.proj):
                epw.write(" proj({})     = {}\n".format(idx+1, pj))
                
            epw.write(" wdata({})    = 'dis_num_iter = {}'\n".format(idx+1, self.epw_inputpara.dis_num_iter))
            epw.write("\n")
            epw.write(" iverbosity  = 2     ! 2 = verbose output for the superconducting part only.\n")
            epw.write(" fsthick  = {}       ! Fermi window thickness [eV]\n".format(self.epw_inputpara.fsthick))
            epw.write(" degaussw = {}       ! smearing in energy-conserving delta functions in [eV]\n".format(self.epw_inputpara.degaussw))
            epw.write(" degaussq = {}       ! smearing for sum over q in the e-ph coupling in [meV]\n".format(self.epw_inputpara.degaussq))
            epw.write(" ephwrite = .false.  ! write ephmatXX, egnv, freq, and ikmap files in prefix.ephmat directory\n")
            epw.write("\n")
            epw.write(" eliashberg     = .true.  ! calculate Eliashberg spectral function\n")
            epw.write(" liso           = .true.  ! solve isotropic ME eqs.\n")
            epw.write(" limag          = .true.  ! solve ME eqs. on imaginary axis\n")
            epw.write(" lpade          = .true.  ! solve ME eqs. on real axis using Pade approximants\n")
            epw.write(" lacon          = .true.  ! analytic continuation of ME eqs. from imaginary to real axis\n")
            epw.write(" nsiter         = 500     ! number of self-consistent iterations when solving ME eqs.\n")
            epw.write(" npade          = {}      ! percentage of Matsubara points used in Pade continuation.\n".format(self.epw_inputpara.npade))
            epw.write(" conv_thr_iaxis = 1.0d-3  ! convergence threshold for solving ME eqs. on imaginary axis\n")
            epw.write(" conv_thr_racon = 1.0d-3  ! convergence threshold for solving ME eqs. on real axis\n")
            epw.write(" wscut          = {}      ! upper limit over Matsubara freq. summation in ME eqs on imag.axis in [eV]\n".format(self.epw_inputpara.wscut))
            epw.write(" muc            = {}      ! effective Coulomb potential used in the ME eqs.\n".format(self.epw_inputpara.muc))
            if  len(self.epw_inputpara.temps) == 2 and \
                float(self.epw_inputpara.temps[0]) < float(self.epw_inputpara.temps[1]) and \
                self.epw_inputpara.nstemp > 2:
                epw.write(" temps      = {} {}   ! number of scf temperature\n".format(self.epw_inputpara.temps[0], self.epw_inputpara.temps[1]))
                epw.write(" nstemp     = {}".format(self.epw_inputpara.nstemp))
            elif len(self.epw_inputpara.temps) > 2 and len(self.epw_inputpara.temps) <= 50:
                epw.write(" temps      = {}      ! number of scf temperature\n".format('  '.join(self.epw_inputpara.temps)))
            elif len(self.epw_inputpara.temps) == 1:
                epw.write(" temps      = {}      ! number of scf temperature\n".format(self.epw_inputpara.temps[0]))
                
                
            #TODO 温度的读取还没写！！！！
            epw.write("\n")
            
            epw.write("\n")
            epw.write(" nk1 = {}\n nk2 = {}\n nk3 = {}\n".format(self.epw_inputpara.nk[0], self.epw_inputpara.nk[1], self.epw_inputpara.nk[2]))
            epw.write("\n")
            epw.write(" nq1 = {}\n nq2 = {}\n nq3 = {}\n".format(self.epw_inputpara.nq[0], self.epw_inputpara.nq[1], self.epw_inputpara.nq[2]))
            epw.write("\n")
            epw.write(" nkf1 = {}\n nkf2 = {}\n nkf3 = {}\n".format(self.epw_inputpara.nkf[0], self.epw_inputpara.nkf[1], self.epw_inputpara.nkf[2]))
            epw.write("\n")
            epw.write(" nqf1 = {}\n nqf2 = {}\n nqf3 = {}\n".format(self.epw_inputpara.nqf[0], self.epw_inputpara.nqf[1], self.epw_inputpara.nqf[2]))
            epw.write("/\n")                                                                

        return inputfilename

    
    def write_epw_aniso_sc_in(self, work_directory:Path):
        inputfilename = "epw_aniso_sc.in"
        epw_phonodata_in = work_directory.joinpath(inputfilename)
        with open(epw_phonodata_in, "w") as epw:
            
            epw.write("&inputepw\n")
            epw.write(" prefix      ='{}',\n".format(self.epw_inputpara.system_name))
            for i, species_name in enumerate(self.epw_inputpara.composition.keys()):
                element      = Element(species_name)
                species_mass = str(element.atomic_mass).strip("amu")
                epw.write(" amass({})    ={},\n".format(i+1, species_mass))          
            epw.write(" outdir      ='./tmp'\n")
            epw.write(" dvscf_dir   ='{}'\n".format(self.epw_inputpara.dvscf_dir))
            epw.write("\n")
            epw.write(" etf_mem     ={}\n".format(self.epw_inputpara.etf_mem))
            epw.write("\n")
            epw.write(" ep_coupling = .false.\n")
            epw.write(" elph        = .false.\n")
            epw.write(" epbwrite    = .false.  ! electron-phonon matrix elements in the coarse Bloch representation and relevant data (dyn matrices)\n")
            epw.write(" epbread     = .false.  ! electron-phonon matrix elements in the coarse Bloch representation and relevant data (dyn matrices)\n")
            epw.write(" epwwrite    = .false.  ! electron-phonon matrix elements in the coarse Wannier representation and relevant data (dyn matrices)\n")
            epw.write(" epwread     = .true.   ! electron-phonon matrix elements in the coarse Wannier representation and relevant data (dyn matrices). ! It is used for a restart calculation and requires kmaps = .true.  prefix.epmatwp, crystal.fmt, vmedata.fmt, wigner.fmt, and epwdata.fmt are all needed for restarting a calculation.\n")
            epw.write("\n")
            epw.write("\n")
            epw.write(" lifc = {}\n".format(self.epw_inputpara.lifc))
            epw.write("\n")
            epw.write(" use_ws      = .true.\n")
            epw.write(" wannierize  = {}\n".format(self.epw_inputpara.wannierize))
            epw.write(" nbndsub     = {}\n".format(self.epw_inputpara.nbndsub))
            if self.epw_inputpara.bands_skipped is not None:
                epw.write(" bands_skipped = 'exclude_bands = {}'\n".format(self.epw_inputpara.bands_skipped))
            epw.write(" num_iter    = {}\n".format(self.epw_inputpara.num_iter))
            epw.write("\n")
            epw.write(" dis_win_max = {} \n".format(self.epw_inputpara.dis_win_max))
            epw.write(" dis_froz_min = {}\n".format(self.epw_inputpara.dis_froz_min))
            epw.write(" dis_froz_max = {}\n".format(self.epw_inputpara.dis_froz_max))
            epw.write("\n")
            for idx, pj in enumerate(self.epw_inputpara.proj):
                epw.write(" proj({})     = {}\n".format(idx+1, pj))
                
            epw.write(" wdata({})    = 'dis_num_iter = {}'\n".format(idx+1, self.epw_inputpara.dis_num_iter))
            epw.write("\n")
            
            epw.write(" iverbosity  = 2     ! 2 = verbose output for the superconducting part only.\n")
            epw.write(" fsthick  = {}       ! Fermi window thickness [eV]\n".format(self.epw_inputpara.fsthick))
            epw.write(" degaussw = {}       ! smearing in energy-conserving delta functions in [eV]\n".format(self.epw_inputpara.degaussw))
            epw.write(" degaussq = {}       ! smearing for sum over q in the e-ph coupling in [meV]\n".format(self.epw_inputpara.degaussq))
            epw.write(" ephwrite = .false.  ! write ephmatXX, egnv, freq, and ikmap files in prefix.ephmat directory\n")
            epw.write("\n")
            
            epw.write(" eliashberg     = .true.  ! calculate Eliashberg spectral function\n")
            epw.write(" laniso         = .true.  ! solve anisotropic ME eqs.\n")
            epw.write(" limag          = .true.  ! solve ME eqs. on imaginary axis\n")
            epw.write(" lpade          = .true.  ! solve ME eqs. on real axis using Pade approximants\n")
            epw.write(" lacon          = .true.  ! analytic continuation of ME eqs. from imaginary to real axis\n")
            epw.write(" nsiter         = 500     ! number of self-consistent iterations when solving ME eqs.\n")
            epw.write(" npade          = {}      ! percentage of Matsubara points used in Pade continuation.\n".format(self.epw_inputpara.npade))
            epw.write(" conv_thr_iaxis = 1.0d-3  ! convergence threshold for solving ME eqs. on imaginary axis\n")
            epw.write(" conv_thr_racon = 1.0d-3  ! convergence threshold for solving ME eqs. on real axis\n")
            epw.write(" wscut          = {}      ! upper limit over Matsubara freq. summation in ME eqs on imag.axis in [eV]\n".format(self.epw_inputpara.wscut))
            epw.write(" muc            = {}      ! effective Coulomb potential used in the ME eqs.\n".format(self.epw_inputpara.muc))
            if  len(self.epw_inputpara.temps) == 2 and \
                float(self.epw_inputpara.temps[0]) < float(self.epw_inputpara.temps[1]) and \
                self.epw_inputpara.nstemp > 2:
                epw.write(" temps      = {} {}   ! number of scf temperature\n".format(self.epw_inputpara.temps[0], self.epw_inputpara.temps[1]))
                epw.write(" nstemp     = {}".format(self.epw_inputpara.nstemp))
            elif len(self.epw_inputpara.temps) > 2 and len(self.epw_inputpara.temps) <= 50:
                epw.write(" temps      = {}      ! number of scf temperature\n".format('  '.join(self.epw_inputpara.temps)))
            elif len(self.epw_inputpara.temps) == 1:
                epw.write(" temps      = {}      ! number of scf temperature\n".format(self.epw_inputpara.temps[0]))
                
                
            #TODO 温度的读取还没写！！！！
            epw.write("\n")
            
            epw.write("\n")
            epw.write(" nk1 = {}\n nk2 = {}\n nk3 = {}\n".format(self.epw_inputpara.nk[0], self.epw_inputpara.nk[1], self.epw_inputpara.nk[2]))
            epw.write("\n")
            epw.write(" nq1 = {}\n nq2 = {}\n nq3 = {}\n".format(self.epw_inputpara.nq[0], self.epw_inputpara.nq[1], self.epw_inputpara.nq[2]))
            epw.write("\n")
            epw.write(" nkf1 = {}\n nkf2 = {}\n nkf3 = {}\n".format(self.epw_inputpara.nkf[0], self.epw_inputpara.nkf[1], self.epw_inputpara.nkf[2]))
            epw.write("\n")
            epw.write(" nqf1 = {}\n nqf2 = {}\n nqf3 = {}\n".format(self.epw_inputpara.nqf[0], self.epw_inputpara.nqf[1], self.epw_inputpara.nqf[2]))
            epw.write("/\n")                                                                

        return inputfilename

