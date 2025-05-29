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
            inputfilename1 = self.write_epw_phono_in(self.epw_inputpara.work_path)
            inputfilename2 = self.write_epw_phonodata_plot_in(self.epw_inputpara.work_path)
            return inputfilename1, inputfilename2
        if mode == "epw_phonodata":
            inputfilename1 = self.write_epw_phonodata_in(self.epw_inputpara.work_path)
            inputfilename2 = self.write_epw_phonodata_plot_in(self.epw_inputpara.work_path)
            return inputfilename
        if mode == "epw_elph":
            inputfilename = self.write_epw_elph_in(self.epw_inputpara.work_path)
            return inputfilename
        if mode == "epw_sc":
            if len(self.epw_inputpara.muc) > 1:
                for mu in self.epw_inputpara.muc:
                    iso_mu_path = self.epw_inputpara.work_path.joinpath("iso_muc_{}".format(mu))
                    aniso_mu_path = self.epw_inputpara.work_path.joinpath("aniso_muc_{}".format(mu))
                    os.makedirs(iso_mu_path, exist_ok=True)
                    os.makedirs(aniso_mu_path, exist_ok=True)
                    inputfilename1 = self.write_epw_iso_sc_in(iso_mu_path, muc=mu)
                    inputfilename2 = self.write_epw_aniso_sc_in(aniso_mu_path, muc=mu)
                return inputfilename1, inputfilename2
            else:
                iso_mu_path = self.epw_inputpara.work_path.joinpath("iso_muc_{}".format(self.epw_inputpara.muc))
                aniso_mu_path = self.epw_inputpara.work_path.joinpath("aniso_muc_{}".format(self.epw_inputpara.muc))
                os.makedirs(iso_mu_path, exist_ok=True)
                os.makedirs(aniso_mu_path, exist_ok=True)
                inputfilename1 = self.write_epw_iso_sc_in(iso_mu_path, muc=self.epw_inputpara.muc)
                inputfilename2 = self.write_epw_aniso_sc_in(aniso_mu_path, muc=self.epw_inputpara.muc)
                return inputfilename1, inputfilename2
        if mode == "epw_prtgkk":
            logger.info("the electron-phonon matrix elements:")
            logger.info(f"the initial electronic states only at k = Gamma(nkf1=nkf2=nkf3=1,  filqf=modified_{self.epw_inputpara.system_name}_band.kpt)")
            prtgkk_path = self.epw_inputpara.work_path.joinpath("prtgkk")
            os.makedirs(prtgkk_path, exist_ok=True)
            inputfilename = self.write_epw_prtgkk_in(prtgkk_path)
            return inputfilename
        if mode == "epw_fermi_nest":
            logger.info("calculate fermi nest:")
            fermi_nest = self.epw_inputpara.work_path.joinpath("fermi_nest")
            os.makedirs(fermi_nest, exist_ok=True)
            inputfilename = self.write_epw_fermi_nest_in(fermi_nest)     
            return inputfilename
        if mode == "epw_linearized_iso":
            if len(self.epw_inputpara.muc) > 1:
                for mu in self.epw_inputpara.muc:
                    iso_mu_path = self.epw_inputpara.work_path.joinpath("iso_muc_{}".format(mu))
                    os.makedirs(iso_mu_path, exist_ok=True)
                    inputfilename1 = self.epw_linearized_iso(iso_mu_path, muc=mu)
                return inputfilename1
            else:
                iso_mu_path = self.epw_inputpara.work_path.joinpath("iso_muc_{}".format(self.epw_inputpara.muc))
                os.makedirs(iso_mu_path, exist_ok=True)
                inputfilename1 = self.epw_linearized_iso(iso_mu_path, muc=self.epw_inputpara.muc)
                return inputfilename1

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
            epw.write(" dvscf_dir   ='{}'      ! directory where .dyn, .dvscf and prefix.phsave/patterns.xx.yy'\n".format(self.epw_inputpara.dvscf_dir))
            epw.write("\n")
            epw.write(" etf_mem     = {}       ! 1 include the whole band 3 only in the window band will be calculated \n".format(self.epw_inputpara.etf_mem))
            epw.write("\n")
            epw.write(" use_ws      = .true.\n")
            epw.write(" wannierize  = .true.   ! If .false., read *.ukk file\n")
            epw.write(" nbndsub     = {}\n".format(self.epw_inputpara.nbndsub))
            if self.epw_inputpara.bands_skipped is not None:
                epw.write(" bands_skipped = 'exclude_bands = {}'\n".format(self.epw_inputpara.bands_skipped))
            epw.write(" num_iter    = {}\n".format(self.epw_inputpara.num_iter))
            epw.write("\n")
            epw.write(" dis_win_max  = {}\n".format(self.epw_inputpara.dis_win_max))
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
            epw.write(" fsthick  = {}       ! Fermi window [eV] : consider only states within Fermi energy +- fsthick\n".format(self.epw_inputpara.fsthick))
            epw.write("\n")
            epw.write(" wdata({})    = 'dis_num_iter = {}'\n".format(idx+5, self.epw_inputpara.dis_num_iter))
            epw.write("\n")
            epw.write(" nk1 = {}\n nk2 = {}\n nk3 = {}\n".format(self.epw_inputpara.nk[0], self.epw_inputpara.nk[1], self.epw_inputpara.nk[2]))
            epw.write("/                           \n")                                                                
        
        return inputfilename
    
    def write_epw_phono_in(self, work_directory:Path):
        inputfilename = "epw_phono.in"
        epw_phonodata_in = work_directory.joinpath(inputfilename)
        with open(epw_phonodata_in, "w") as epw:
            
            epw.write("&inputepw\n")
            epw.write(" prefix      ='{}',\n".format(self.epw_inputpara.system_name))
            for i, species_name in enumerate(self.epw_inputpara.composition.keys()):
                element      = Element(species_name)
                species_mass = str(element.atomic_mass).strip("amu")
                epw.write(" amass({})    ={},\n".format(i+1, species_mass))          
            epw.write(" outdir      ='./tmp'\n")
            epw.write(" dvscf_dir   ='{}'      ! directory where .dyn, .dvscf and prefix.phsave/patterns.xx.yy'\n".format(self.epw_inputpara.dvscf_dir))
            epw.write("\n")
            epw.write(" etf_mem     = {}       ! 1 include the whole band 3 only in the window band will be calculated \n".format(self.epw_inputpara.etf_mem))
            epw.write("\n")
            epw.write(" elph        = .true.   ! If .true. calculate e-ph coefficients.\n")
            epw.write(" epbwrite    = .true.   ! electron-phonon matrix elements in the coarse Bloch representation and relevant data (dyn matrices)\n")
            epw.write(" epbread     = .false.  ! electron-phonon matrix elements in the coarse Bloch representation and relevant data (dyn matrices)\n")
            epw.write(" epwwrite    = .true.   ! electron-phonon matrix elements in the coarse Wannier representation and relevant data (dyn matrices)\n")
            epw.write(" epwread     = .false.  ! electron-phonon matrix elements in the coarse Wannier representation and relevant data (dyn matrices). ! It is used for a restart calculation and requires kmaps = .true.  prefix.epmatwp, crystal.fmt, vmedata.fmt, wigner.fmt, and epwdata.fmt are all needed for restarting a calculation.\n")
            
            epw.write("\n")
            epw.write(" lifc = {}\n".format(self.epw_inputpara.lifc))
            epw.write("\n")
            
            epw.write(" use_ws      = .true.\n")
            epw.write(" wannierize  = {}       ! If .false., read *.ukk file\n".format(self.epw_inputpara.wannierize))
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
            epw.write(" fsthick  = {}       ! Fermi window [eV] : consider only states within Fermi energy +- fsthick\n".format(self.epw_inputpara.fsthick))
            epw.write("\n")
            epw.write(" band_plot    = .true.   \n")
            epw.write(" filkf        = '{}'     \n".format(self.epw_inputpara.filkf))
            epw.write(" filqf        = '{}'     \n".format(self.epw_inputpara.filqf))
            epw.write(" asr_typ      = '{}'     \n".format(self.epw_inputpara.asr_typ))

            epw.write("\n")
            epw.write(" nk1 = {}\n nk2 = {}\n nk3 = {}\n".format(self.epw_inputpara.nk[0], self.epw_inputpara.nk[1], self.epw_inputpara.nk[2]))
            epw.write("\n")
            epw.write(" nq1 = {}\n nq2 = {}\n nq3 = {}\n".format(self.epw_inputpara.nq[0], self.epw_inputpara.nq[1], self.epw_inputpara.nq[2]))
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
            epw.write(" dvscf_dir   ='{}'      ! directory where .dyn, .dvscf and prefix.phsave/patterns.xx.yy'\n".format(self.epw_inputpara.dvscf_dir))
            epw.write("\n")
            epw.write(" etf_mem     = {}       ! 1 include the whole band 3 only in the window band will be calculated \n".format(self.epw_inputpara.etf_mem))
            epw.write("\n")
            epw.write(" elph        = .true.   ! calculate e-ph coefficients.   ephwrite requires ep_coupling=.TRUE., elph=.TRUE.\n")
            epw.write(" epbwrite    = .false.  ! electron-phonon matrix elements in the coarse Bloch representation and relevant data (dyn matrices)\n")
            epw.write(" epbread     = .false.  ! electron-phonon matrix elements in the coarse Bloch representation and relevant data (dyn matrices)\n")
            epw.write(" epwwrite    = .false.  ! write electron-phonon matrix elements in the coarse Wannier representation and relevant data (dyn matrices)\n")
            epw.write(" epwread     = .true.   ! read e-ph matrices from 'prefix.epmatwp' file. electron-phonon matrix elements in the coarse Wannier representation and relevant data (dyn matrices). \n")
            epw.write("                        ! It is used for a restart calculation, and doesn't set kmaps = .true. kmaps may appear in EPW website, it's out of date.  prefix.epmatwp, crystal.fmt, vmedata.fmt, wigner.fmt, and epwdata.fmt are all needed for restarting a calculation.\n")
            
            epw.write("\n")
            epw.write(" lifc = {}\n".format(self.epw_inputpara.lifc))
            epw.write("\n")
            
            epw.write(" use_ws      = .true.\n")
            epw.write(" wannierize  = {}       ! If .false., read *.ukk file\n".format(self.epw_inputpara.wannierize))
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
            epw.write(" fsthick  = {}       ! Fermi window [eV] : consider only states within Fermi energy +- fsthick\n".format(self.epw_inputpara.fsthick))
            epw.write("\n")
            epw.write(" band_plot    = .true.   \n")
            epw.write(" filkf        = '{}'     \n".format(self.epw_inputpara.filkf))
            epw.write(" filqf        = '{}'     \n".format(self.epw_inputpara.filqf))
            epw.write(" asr_typ      = '{}'     \n".format(self.epw_inputpara.asr_typ))

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
            epw.write("0 4000\n")
            epw.write("freq.dat\n")
            epw.write("freq.dat.ps\n")
            epw.write("0\n")
            epw.write("0 50\n")
            epw.write("\n\n\n")
            # epw.write("# Second line is min and max value in your spectrum\n")
            # epw.write("# Meanwhile, another file named band.eig is eletron spectrum.\n")
            # epw.write("# band.eig also can be plotted by plotband.x, detailed in Tue.6.Lafuente.pdf.\n") 
        return inputfilename

    def write_epw_elph_in(self, work_directory:Path):
        inputfilename = "epw_elph.in"
        epw_elph_in = work_directory.joinpath(inputfilename)
        with open(epw_elph_in, "w") as epw:
            
            epw.write("&inputepw\n")
            epw.write(" prefix      ='{}',\n".format(self.epw_inputpara.system_name))
            for i, species_name in enumerate(self.epw_inputpara.composition.keys()):
                element      = Element(species_name)
                species_mass = str(element.atomic_mass).strip("amu")
                epw.write(" amass({})    ={},\n".format(i+1, species_mass))          
            epw.write(" outdir      ='./tmp'\n")
            epw.write(" dvscf_dir   ='{}'      ! directory where .dyn, .dvscf and prefix.phsave/patterns.xx.yy'\n".format(self.epw_inputpara.dvscf_dir))
            epw.write("\n")
            epw.write(" etf_mem     = {}       ! 1 include the whole band 3 only in the window band will be calculated \n".format(self.epw_inputpara.etf_mem))
            epw.write("\n")
            epw.write(" ep_coupling = {}       ! run e-ph coupling calculation. ephwrite requires ep_coupling=.TRUE., elph=.TRUE.\n".format(self.epw_inputpara.ep_coupling))
            epw.write(" elph        = {}       ! calculate e-ph coefficients.   ephwrite requires ep_coupling=.TRUE., elph=.TRUE.\n".format(self.epw_inputpara.elph))
            epw.write(" epbwrite    = {}       ! electron-phonon matrix elements in the coarse Bloch representation and relevant data (dyn matrices)\n".format(self.epw_inputpara.epbwrite))
            epw.write(" epbread     = {}       ! electron-phonon matrix elements in the coarse Bloch representation and relevant data (dyn matrices)\n".format(self.epw_inputpara.epbread))
            epw.write(" epwwrite    = {}       ! write electron-phonon matrix elements in the coarse Wannier representation and relevant data (dyn matrices)\n".format(self.epw_inputpara.epwwrite))
            epw.write(" epwread     = {}       ! read e-ph matrices from 'prefix.epmatwp' file. electron-phonon matrix elements in the coarse Wannier representation and relevant data (dyn matrices). \n".format(self.epw_inputpara.epwread))
            epw.write("                        ! It is used for a restart calculation, and doesn't set kmaps = .true. kmaps may appear in EPW website, it's out of date.  prefix.epmatwp, crystal.fmt, vmedata.fmt, wigner.fmt, and epwdata.fmt are all needed for restarting a calculation.\n")
            epw.write(" ephwrite    = {}       ! write ephmatXX, egnv, freq, and ikmap files in prefix.ephmat directory\n".format(self.epw_inputpara.ephwrite))
            epw.write("                        ! ephmatXX files (one per CPU) containing the electron-phonon matrix elements within the Fermi window (fsthick) on the dense k and q grids\n")
            epw.write("                        ! freq file containing the phonon frequencies on the dense q grid\n")
            epw.write("                        ! egnv file containing the eigenvalues within the Fermi window on the dense k grid\n")
            epw.write("                        ! ikmap file containing the index of the k-points on the dense (irreducible) grid within the Fermi window\n")

            epw.write("\n")
            epw.write(" lifc = {}\n".format(self.epw_inputpara.lifc))
            epw.write("\n")
            
            epw.write(" use_ws      = .true.\n")
            epw.write(" wannierize  = {}       ! If .false., read *.ukk file\n".format(self.epw_inputpara.wannierize))
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
            epw.write(" fsthick  = {}       ! Fermi window [eV] : consider only states within Fermi energy +- fsthick\n".format(self.epw_inputpara.fsthick))
            epw.write(" degaussw = {}       ! smearing in energy-conserving delta functions in [eV]\n".format(self.epw_inputpara.degaussw))
            epw.write(" degaussq = {}       ! smearing for sum over q in the e-ph coupling in [meV]\n".format(self.epw_inputpara.degaussq))
            epw.write(" delta_qsmear = {} \n".format(self.epw_inputpara.delta_qsmear))
            epw.write(" nqsmear = {} \n".format(self.epw_inputpara.nqsmear))    

            epw.write("\n")
            epw.write(" nk1 = {}\n nk2 = {}\n nk3 = {}\n".format(self.epw_inputpara.nk[0], self.epw_inputpara.nk[1], self.epw_inputpara.nk[2]))
            epw.write("\n")
            epw.write(" nq1 = {}\n nq2 = {}\n nq3 = {}\n".format(self.epw_inputpara.nq[0], self.epw_inputpara.nq[1], self.epw_inputpara.nq[2]))
            epw.write("\n")
            epw.write(" mp_mesh_k = .true.\n")
            epw.write(" nkf1 = {}\n nkf2 = {}\n nkf3 = {}\n".format(self.epw_inputpara.nkf[0], self.epw_inputpara.nkf[1], self.epw_inputpara.nkf[2]))
            epw.write("\n")
            epw.write(" nqf1 = {}\n nqf2 = {}\n nqf3 = {}\n".format(self.epw_inputpara.nqf[0], self.epw_inputpara.nqf[1], self.epw_inputpara.nqf[2]))
            epw.write("/                           \n")                                                                

        return inputfilename

    def write_epw_iso_sc_in(self, work_directory:Path, muc=None):
        inputfilename = "epw_iso_sc.in"
        epw_iso_sc_in = work_directory.joinpath(inputfilename)
        with open(epw_iso_sc_in, "w") as epw:
            
            epw.write("&inputepw\n")
            epw.write(" prefix      ='{}',\n".format(self.epw_inputpara.system_name))
            for i, species_name in enumerate(self.epw_inputpara.composition.keys()):
                element      = Element(species_name)
                species_mass = str(element.atomic_mass).strip("amu")
                epw.write(" amass({})    ={},\n".format(i+1, species_mass))          
            epw.write(" outdir      ='./tmp'\n")
            epw.write(" dvscf_dir   ='{}'      ! directory where .dyn, .dvscf and prefix.phsave/patterns.xx.yy'\n".format(self.epw_inputpara.dvscf_dir))
            epw.write("\n")
            epw.write(" etf_mem     = {}       ! 1 include the whole band 3 only in the window band will be calculated \n".format(self.epw_inputpara.etf_mem))
            epw.write("\n")
            epw.write(" ep_coupling = {}       ! run e-ph coupling calculation. ephwrite requires ep_coupling=.TRUE., elph=.TRUE.\n".format(self.epw_inputpara.ep_coupling))
            epw.write(" elph        = {}       ! calculate e-ph coefficients.   ephwrite requires ep_coupling=.TRUE., elph=.TRUE.\n".format(self.epw_inputpara.elph))
            epw.write(" epbwrite    = {}       ! electron-phonon matrix elements in the coarse Bloch representation and relevant data (dyn matrices)\n".format(self.epw_inputpara.epbwrite))
            epw.write(" epbread     = {}       ! electron-phonon matrix elements in the coarse Bloch representation and relevant data (dyn matrices)\n".format(self.epw_inputpara.epbread))
            epw.write(" epwwrite    = {}       ! write electron-phonon matrix elements in the coarse Wannier representation and relevant data (dyn matrices)\n".format(self.epw_inputpara.epwwrite))
            epw.write(" epwread     = {}       ! read e-ph matrices from 'prefix.epmatwp' file. electron-phonon matrix elements in the coarse Wannier representation and relevant data (dyn matrices). \n".format(self.epw_inputpara.epwread))
            epw.write("                        ! It is used for a restart calculation, and doesn't set kmaps = .true. kmaps may appear in EPW website, it's out of date.  prefix.epmatwp, crystal.fmt, vmedata.fmt, wigner.fmt, and epwdata.fmt are all needed for restarting a calculation.\n")
            epw.write(" ephwrite    = {}       ! write ephmatXX, egnv, freq, and ikmap files in prefix.ephmat directory\n".format(self.epw_inputpara.ephwrite))
            epw.write("                        ! ephmatXX files (one per CPU) containing the electron-phonon matrix elements within the Fermi window (fsthick) on the dense k and q grids\n")
            epw.write("                        ! freq file containing the phonon frequencies on the dense q grid\n")
            epw.write("                        ! egnv file containing the eigenvalues within the Fermi window on the dense k grid\n")
            epw.write("                        ! ikmap file containing the index of the k-points on the dense (irreducible) grid within the Fermi window\n")

            epw.write("\n")
            epw.write(" lifc = {}\n".format(self.epw_inputpara.lifc))
            epw.write("\n")
            
            epw.write(" use_ws      = .true.\n")
            epw.write(" wannierize  = {}       ! If .false., read *.ukk file\n".format(self.epw_inputpara.wannierize))
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
            epw.write("\n")
            epw.write(" fsthick  = {}        ! Fermi window [eV] : consider only states within Fermi energy +- fsthick\n".format(self.epw_inputpara.fsthick))
            epw.write(" degaussw = {}        ! smearing in energy-conserving delta functions in [eV]\n".format(self.epw_inputpara.degaussw))
            epw.write(" degaussq = {}        ! smearing for sum over q in the e-ph coupling in [meV]\n".format(self.epw_inputpara.degaussq))
            epw.write(" delta_qsmear = {} \n".format(self.epw_inputpara.delta_qsmear))
            epw.write(" nqsmear = {} \n".format(self.epw_inputpara.nqsmear))    
            epw.write(" selecqread = .false. ! If .true. then restart from the selecq.fmt file.\n")
            epw.write("                      ! selecq.fmt only needed if selecqread=.true. otherwise it will be re-created")
            epw.write("\n")
            
            epw.write(" eliashberg     = .true.  ! calculate Eliashberg spectral function\n")
            epw.write(" liso           = .true.  ! solve isotropic ME eqs.\n")
            epw.write(" limag          = .true.  ! solve ME eqs. on imaginary axis\n")
            epw.write(" lpade          = .true.  ! solve ME eqs. on real axis using Pade approximants\n")
            epw.write(" lacon          = {}      ! analytic continuation of ME eqs. from imaginary to real axis. it takes very long time, if you need it, you can turn it on!\n".format(self.epw_inputpara.lacon))
            epw.write(" nsiter         = 500     ! number of self-consistent iterations when solving ME eqs.\n")
            epw.write(" npade          = {}      ! percentage of Matsubara points used in Pade continuation.\n".format(self.epw_inputpara.npade))
            epw.write(" conv_thr_iaxis = 1.0d-3  ! convergence threshold for solving ME eqs. on imaginary axis\n")
            epw.write(" conv_thr_racon = 1.0d-3  ! convergence threshold for solving ME eqs. on real axis\n")
            epw.write(" wscut          = {}      ! upper limit over Matsubara freq. summation in ME eqs on imag.axis in [eV]\n".format(self.epw_inputpara.wscut))
            epw.write(" muc            = {}      ! effective Coulomb potential used in the ME eqs.\n".format(muc))
            if  len(self.epw_inputpara.temps) == 2 and \
                float(self.epw_inputpara.temps[0]) < float(self.epw_inputpara.temps[1]) and \
                self.epw_inputpara.nstemp > 2:
                epw.write(" temps      = {} {}   \n".format(self.epw_inputpara.temps[0], self.epw_inputpara.temps[1]))
                epw.write(" nstemp     = {}      ! number of scf temperature".format(self.epw_inputpara.nstemp))
            elif len(self.epw_inputpara.temps) > 2 and len(self.epw_inputpara.temps) <= 50:
                epw.write(" temps      = {}      ! many values of scf temperature\n".format('  '.join(self.epw_inputpara.temps)))
            elif len(self.epw_inputpara.temps) == 1:
                epw.write(" temps      = {}      ! only one value scf temperature\n".format(self.epw_inputpara.temps[0]))
                
                
            epw.write("\n")
            
            epw.write("\n")
            epw.write(" nk1 = {}\n nk2 = {}\n nk3 = {}\n".format(self.epw_inputpara.nk[0], self.epw_inputpara.nk[1], self.epw_inputpara.nk[2]))
            epw.write("\n")
            epw.write(" nq1 = {}\n nq2 = {}\n nq3 = {}\n".format(self.epw_inputpara.nq[0], self.epw_inputpara.nq[1], self.epw_inputpara.nq[2]))
            epw.write("\n")
            epw.write(" mp_mesh_k = .true.\n")
            epw.write(" nkf1 = {}\n nkf2 = {}\n nkf3 = {}\n".format(self.epw_inputpara.nkf[0], self.epw_inputpara.nkf[1], self.epw_inputpara.nkf[2]))
            epw.write("\n")
            epw.write(" nqf1 = {}\n nqf2 = {}\n nqf3 = {}\n".format(self.epw_inputpara.nqf[0], self.epw_inputpara.nqf[1], self.epw_inputpara.nqf[2]))
            epw.write("/\n")                                                                

        return inputfilename

    def epw_linearized_iso(self, work_directory:Path, muc=None):
        inputfilename = "epw_linearized_iso.in"
        epw_linearized_iso_in = work_directory.joinpath(inputfilename)
        with open(epw_linearized_iso_in, "w") as epw:

            epw.write("&inputepw\n")
            epw.write(" prefix      ='{}',\n".format(self.epw_inputpara.system_name))
            for i, species_name in enumerate(self.epw_inputpara.composition.keys()):
                element      = Element(species_name)
                species_mass = str(element.atomic_mass).strip("amu")
                epw.write(" amass({})    ={},\n".format(i+1, species_mass))          
            epw.write(" outdir      ='./tmp'\n")
            epw.write(" dvscf_dir   ='{}'           ! directory where .dyn, .dvscf and prefix.phsave/patterns.xx.yy'\n".format(self.epw_inputpara.dvscf_dir))
            epw.write("\n")
            epw.write(" etf_mem     = 1             ! 1 include the whole band 3 only in the window band will be calculated \n")
            epw.write("\n")
            epw.write(" ep_coupling = .false.       ! run e-ph coupling calculation. ephwrite requires ep_coupling=.TRUE., elph=.TRUE.\n")
            epw.write(" elph        = .false.       ! calculate e-ph coefficients.   ephwrite requires ep_coupling=.TRUE., elph=.TRUE.\n")
            epw.write(" epbwrite    = .false.       ! electron-phonon matrix elements in the coarse Bloch representation and relevant data (dyn matrices)\n")
            epw.write(" epbread     = .false.       ! electron-phonon matrix elements in the coarse Bloch representation and relevant data (dyn matrices)\n")
            epw.write(" epwwrite    = .false.       ! write electron-phonon matrix elements in the coarse Wannier representation and relevant data (dyn matrices)\n")
            epw.write(" epwread     = .true.        ! read e-ph matrices from 'prefix.epmatwp' file. electron-phonon matrix elements in the coarse Wannier representation and relevant data (dyn matrices). \n")
            epw.write("                             ! It is used for a restart calculation, and doesn't set kmaps = .true. kmaps may appear in EPW website, it's out of date.  prefix.epmatwp, crystal.fmt, vmedata.fmt, wigner.fmt, and epwdata.fmt are all needed for restarting a calculation.\n")
            epw.write(" ephwrite    = .false.       ! write ephmatXX, egnv, freq, and ikmap files in prefix.ephmat directory\n")
            epw.write("\n")
            
            epw.write("\n")
            epw.write(" lifc = {}\n".format(self.epw_inputpara.lifc))
            epw.write("\n")
            
            epw.write(" use_ws      = .true.\n")
            epw.write(" wannierize  = {}       ! If .false., read *.ukk file\n".format(self.epw_inputpara.wannierize))
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
            epw.write(" iverbosity  = 2      ! 2 = verbose output for the superconducting part only.\n")
            epw.write("\n")
            epw.write(" fsthick  = {}        ! Fermi window [eV] : consider only states within Fermi energy +- fsthick\n".format(self.epw_inputpara.fsthick))
            epw.write(" degaussw = {}        ! smearing in energy-conserving delta functions in [eV]\n".format(self.epw_inputpara.degaussw))
            epw.write(" degaussq = {}        ! smearing for sum over q in the e-ph coupling in [meV]\n".format(self.epw_inputpara.degaussq))
            epw.write(" delta_qsmear = {}    \n".format(self.epw_inputpara.delta_qsmear))
            epw.write(" nqsmear = {}         \n".format(self.epw_inputpara.nqsmear))    
            epw.write(" selecqread = .false. ! If .true. then restart from the selecq.fmt file.\n")
            epw.write("                      ! selecq.fmt only needed if selecqread=.true. otherwise it will be re-created")
            epw.write("\n")
            
            epw.write(" eliashberg       = .true.  ! calculate Eliashberg spectral function\n")
            epw.write(" liso             = .true.  ! solve isotropic ME eqs.\n")
            epw.write(" limag            = .true.  ! solve ME eqs. on imaginary axis\n")
            epw.write(" lpade            = .false. ! solve ME eqs. on real axis using Pade approximants\n")
            epw.write(" lacon            = .false. ! analytic continuation of ME eqs. from imaginary to real axis. it takes very long time, if you need it, you can turn it on!\n".format(self.epw_inputpara.lacon))
            epw.write(" nsiter           = 500     ! number of self-consistent iterations when solving ME eqs.\n")
            epw.write(" npade            = {}      ! percentage of Matsubara points used in Pade continuation.\n".format(self.epw_inputpara.npade))
            epw.write(" conv_thr_iaxis   = 1.0d-3  ! convergence threshold for solving ME eqs. on imaginary axis\n")
            epw.write(" conv_thr_racon   = 1.0d-3  ! convergence threshold for solving ME eqs. on real axis\n")
            epw.write(" fila2f           = {}.a2f  ! read the a2f file\n".format(self.epw_inputpara.system_name))
            epw.write(" tc_linear        = .true.  ! solve linearized ME eqn. for Tc\n")
            epw.write(" tc_linear_solver ='power'  ! algorithm to solve Tc eigenvalue problem: 'power' OR 'lapack'\n")
            epw.write(" wscut            = {}      ! upper limit over Matsubara freq. summation in ME eqs on imag.axis in [eV]\n".format(self.epw_inputpara.wscut))
            epw.write(" muc              = {}      ! effective Coulomb potential used in the ME eqs.\n".format(muc))
            epw.write("\n")

            if  len(self.epw_inputpara.temps) == 2 and \
                float(self.epw_inputpara.temps[0]) < float(self.epw_inputpara.temps[1]) and \
                self.epw_inputpara.nstemp > 2:
                epw.write(" temps      = {} {}   ! number of scf temperature\n".format(self.epw_inputpara.temps[0], self.epw_inputpara.temps[1]))
                epw.write(" nstemp     = {}".format(self.epw_inputpara.nstemp))
            elif len(self.epw_inputpara.temps) > 2 and len(self.epw_inputpara.temps) <= 50:
                epw.write(" temps      = {}      ! number of scf temperature\n".format('  '.join(self.epw_inputpara.temps)))
            elif len(self.epw_inputpara.temps) == 1:
                epw.write(" temps      = {}      ! number of scf temperature\n".format(self.epw_inputpara.temps[0]))

            epw.write("\n")
            epw.write("\n")
            epw.write(" nk1 = {}\n nk2 = {}\n nk3 = {}\n".format(self.epw_inputpara.nk[0], self.epw_inputpara.nk[1], self.epw_inputpara.nk[2]))
            epw.write("\n")
            epw.write(" nq1 = {}\n nq2 = {}\n nq3 = {}\n".format(self.epw_inputpara.nq[0], self.epw_inputpara.nq[1], self.epw_inputpara.nq[2]))
            epw.write("\n")
            epw.write(" mp_mesh_k = .true.\n")
            epw.write(" nkf1 = {}\n nkf2 = {}\n nkf3 = {}\n".format(self.epw_inputpara.nkf[0], self.epw_inputpara.nkf[1], self.epw_inputpara.nkf[2]))
            epw.write("\n")
            epw.write(" nqf1 = {}\n nqf2 = {}\n nqf3 = {}\n".format(self.epw_inputpara.nqf[0], self.epw_inputpara.nqf[1], self.epw_inputpara.nqf[2]))
            epw.write("/\n")                                                                

        return inputfilename

    def write_epw_aniso_sc_in(self, work_directory:Path, muc=None):
        inputfilename = "epw_aniso_sc.in"
        epw_aniso_sc_in = work_directory.joinpath(inputfilename)
        with open(epw_aniso_sc_in, "w") as epw:
            
            epw.write("&inputepw\n")
            epw.write(" prefix      ='{}',\n".format(self.epw_inputpara.system_name))
            for i, species_name in enumerate(self.epw_inputpara.composition.keys()):
                element      = Element(species_name)
                species_mass = str(element.atomic_mass).strip("amu")
                epw.write(" amass({})    ={},\n".format(i+1, species_mass))          
            epw.write(" outdir      ='./tmp'\n")
            epw.write(" dvscf_dir   ='{}'      ! directory where .dyn, .dvscf and prefix.phsave/patterns.xx.yy'\n".format(self.epw_inputpara.dvscf_dir))
            epw.write("\n")
            epw.write(" etf_mem     = {}       ! 1 include the whole band 3 only in the window band will be calculated \n".format(self.epw_inputpara.etf_mem))
            epw.write("\n")
            epw.write(" ep_coupling = {}       ! run e-ph coupling calculation. ephwrite requires ep_coupling=.TRUE., elph=.TRUE.\n".format(self.epw_inputpara.ep_coupling))
            epw.write(" elph        = {}       ! calculate e-ph coefficients.   ephwrite requires ep_coupling=.TRUE., elph=.TRUE.\n".format(self.epw_inputpara.elph))
            epw.write(" epbwrite    = {}       ! electron-phonon matrix elements in the coarse Bloch representation and relevant data (dyn matrices)\n".format(self.epw_inputpara.epbwrite))
            epw.write(" epbread     = {}       ! electron-phonon matrix elements in the coarse Bloch representation and relevant data (dyn matrices)\n".format(self.epw_inputpara.epbread))
            epw.write(" epwwrite    = {}       ! write electron-phonon matrix elements in the coarse Wannier representation and relevant data (dyn matrices)\n".format(self.epw_inputpara.epwwrite))
            epw.write(" epwread     = {}       ! read e-ph matrices from 'prefix.epmatwp' file. electron-phonon matrix elements in the coarse Wannier representation and relevant data (dyn matrices). \n".format(self.epw_inputpara.epwread))
            epw.write("                        ! It is used for a restart calculation, and doesn't set kmaps = .true. kmaps may appear in EPW website, it's out of date.  prefix.epmatwp, crystal.fmt, vmedata.fmt, wigner.fmt, and epwdata.fmt are all needed for restarting a calculation.\n")
            epw.write(" ephwrite    = {}       ! write ephmatXX, egnv, freq, and ikmap files in prefix.ephmat directory\n".format(self.epw_inputpara.ephwrite))
            epw.write("                        ! ephmatXX files (one per CPU) containing the electron-phonon matrix elements within the Fermi window (fsthick) on the dense k and q grids\n")
            epw.write("                        ! freq file containing the phonon frequencies on the dense q grid\n")
            epw.write("                        ! egnv file containing the eigenvalues within the Fermi window on the dense k grid\n")
            epw.write("                        ! ikmap file containing the index of the k-points on the dense (irreducible) grid within the Fermi window\n")

            epw.write("\n")
            epw.write(" lifc = {}\n".format(self.epw_inputpara.lifc))
            epw.write("\n")
            
            epw.write(" use_ws      = .true.\n")
            epw.write(" wannierize  = {}       ! If .false., read *.ukk file\n".format(self.epw_inputpara.wannierize))
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
            
            epw.write(" iverbosity  = 2      ! 2 = verbose output for the superconducting part only.\n")
            epw.write(" fsthick  = {}        ! Fermi window [eV] : consider only states within Fermi energy +- fsthick\n".format(self.epw_inputpara.fsthick))
            epw.write(" degaussw = {}        ! smearing in energy-conserving delta functions in [eV]\n".format(self.epw_inputpara.degaussw))
            epw.write(" degaussq = {}        ! smearing for sum over q in the e-ph coupling in [meV]\n".format(self.epw_inputpara.degaussq))
            epw.write(" delta_qsmear = {} \n".format(self.epw_inputpara.delta_qsmear))
            epw.write(" nqsmear = {} \n".format(self.epw_inputpara.nqsmear))    
            epw.write(" selecqread = .false. ! If .true. then restart from the selecq.fmt file.\n")
            epw.write("                      ! selecq.fmt only needed if selecqread=.true. otherwise it will be re-created")
            epw.write("\n")
            
            epw.write(" eliashberg     = .true.  ! calculate Eliashberg spectral function\n")
            epw.write(" laniso         = .true.  ! solve anisotropic ME eqs.\n")
            epw.write(" limag          = .true.  ! solve ME eqs. on imaginary axis\n")
            epw.write(" lpade          = .true.  ! solve ME eqs. on real axis using Pade approximants\n")
            epw.write(" lacon          = {}      ! analytic continuation of ME eqs. from imaginary to real axis. it takes very long time, if you need it, you can turn it on!\n".format(self.epw_inputpara.lacon))
            epw.write(" nsiter         = 500     ! number of self-consistent iterations when solving ME eqs.\n")
            epw.write(" npade          = {}      ! percentage of Matsubara points used in Pade continuation.\n".format(self.epw_inputpara.npade))
            epw.write(" conv_thr_iaxis = 1.0d-3  ! convergence threshold for solving ME eqs. on imaginary axis\n")
            epw.write(" conv_thr_racon = 1.0d-3  ! convergence threshold for solving ME eqs. on real axis\n")
            epw.write(" wscut          = {}      ! upper limit over Matsubara freq. summation in ME eqs on imag.axis in [eV]\n".format(self.epw_inputpara.wscut))
            epw.write(" muc            = {}      ! effective Coulomb potential used in the ME eqs.\n".format(muc))
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
            epw.write(" mp_mesh_k = .true.\n")
            epw.write(" nkf1 = {}\n nkf2 = {}\n nkf3 = {}\n".format(self.epw_inputpara.nkf[0], self.epw_inputpara.nkf[1], self.epw_inputpara.nkf[2]))
            epw.write("\n")
            epw.write(" nqf1 = {}\n nqf2 = {}\n nqf3 = {}\n".format(self.epw_inputpara.nqf[0], self.epw_inputpara.nqf[1], self.epw_inputpara.nqf[2]))
            epw.write("/\n")                                                                

        return inputfilename

    def write_epw_prtgkk_in(self, work_directory:Path):
        inputfilename = "epw_prtgkk.in"
        epw_prtgkk_in = work_directory.joinpath(inputfilename)
        with open(epw_prtgkk_in, "w") as epw:

            epw.write("&inputepw\n")
            epw.write(" prefix      ='{}',\n".format(self.epw_inputpara.system_name))
            for i, species_name in enumerate(self.epw_inputpara.composition.keys()):
                element      = Element(species_name)
                species_mass = str(element.atomic_mass).strip("amu")
                epw.write(" amass({})    ={},\n".format(i+1, species_mass))          
            epw.write(" outdir      ='./tmp'\n")
            epw.write(" dvscf_dir   ='{}'      ! directory where .dyn, .dvscf and prefix.phsave/patterns.xx.yy'\n".format(self.epw_inputpara.dvscf_dir))
            epw.write("\n")
            epw.write(" etf_mem     = {}       ! 1 include the whole band 3 only in the window band will be calculated \n".format(self.epw_inputpara.etf_mem))
            epw.write("\n")
            epw.write(" elph        = .true.       ! calculate e-ph coefficients.   ephwrite requires ep_coupling=.TRUE., elph=.TRUE.\n")
            epw.write(" epwwrite    = .false.      ! write electron-phonon matrix elements in the coarse Wannier representation and relevant data (dyn matrices)\n")
            epw.write(" epwread     = .true.       ! read e-ph matrices from 'prefix.epmatwp' file. electron-phonon matrix elements in the coarse Wannier representation and relevant data (dyn matrices). \n")
            epw.write("\n")
            
            epw.write(" prtgkk      = .true.     \n")
            
            epw.write(" use_ws      = .true.     \n")
            epw.write(" wannierize  = .false.      ! If .false., read *.ukk file\n".format(self.epw_inputpara.wannierize))
            epw.write("\n")  
            epw.write(" iverbosity  = 2     ! 2 = verbose output for the superconducting part only.\n")
            epw.write("\n")
            epw.write(" nk1 = {}\n nk2 = {}\n nk3 = {}\n".format(self.epw_inputpara.nk[0], self.epw_inputpara.nk[1], self.epw_inputpara.nk[2]))
            epw.write("\n")
            epw.write(" nq1 = {}\n nq2 = {}\n nq3 = {}\n".format(self.epw_inputpara.nq[0], self.epw_inputpara.nq[1], self.epw_inputpara.nq[2]))
            epw.write("\n")

            epw.write(" nkf1 = {}\n nkf2 = {}\n nkf3 = {}\n".format(self.epw_inputpara.nkf[0], self.epw_inputpara.nkf[1], self.epw_inputpara.nkf[2]))
            epw.write("\n")
            epw.write("! nqf1 = {}\n nqf2 = {}\n nqf3 = {}\n".format(self.epw_inputpara.nqf[0], self.epw_inputpara.nqf[1], self.epw_inputpara.nqf[2]))
            epw.write(" filqf        = '{}'     \n".format(self.epw_inputpara.filqf))
            epw.write("/                           \n")                                                                

        return inputfilename

    def write_epw_fermi_nest_in(self, work_directory:Path):
        inputfilename = "epw_fermi_nest.in"
        epw_fermi_nest_in = work_directory.joinpath(inputfilename)
        with open(epw_fermi_nest_in, "w") as epw:

            epw.write("&inputepw\n")
            epw.write(" prefix      ='{}',\n".format(self.epw_inputpara.system_name))
            for i, species_name in enumerate(self.epw_inputpara.composition.keys()):
                element      = Element(species_name)
                species_mass = str(element.atomic_mass).strip("amu")
                epw.write(" amass({})    ={},\n".format(i+1, species_mass))          
            epw.write(" outdir      ='./tmp'\n")
            epw.write(" dvscf_dir   ='{}'      ! directory where .dyn, .dvscf and prefix.phsave/patterns.xx.yy'\n".format(self.epw_inputpara.dvscf_dir))
            epw.write("\n")
            epw.write(" etf_mem     = {}       ! 1 include the whole band 3 only in the window band will be calculated \n".format(self.epw_inputpara.etf_mem))
            epw.write("\n")
            epw.write(" elph        = .true.       ! calculate e-ph coefficients.   ephwrite requires ep_coupling=.TRUE., elph=.TRUE.\n")
            epw.write(" epwwrite    = .false.      ! write electron-phonon matrix elements in the coarse Wannier representation and relevant data (dyn matrices)\n")
            epw.write(" epwread     = .true.       ! read e-ph matrices from 'prefix.epmatwp' file. electron-phonon matrix elements in the coarse Wannier representation and relevant data (dyn matrices). \n")
            epw.write("\n")
            
            epw.write(" wannierize  = {}       ! If .false., read *.ukk file\n".format(self.epw_inputpara.wannierize))
            epw.write(" nbndsub     = {}\n".format(self.epw_inputpara.nbndsub))
            if self.epw_inputpara.bands_skipped is not None:
                epw.write(" bands_skipped = 'exclude_bands = {}'\n".format(self.epw_inputpara.bands_skipped))
            epw.write(" num_iter    = {}\n".format(self.epw_inputpara.num_iter))
            epw.write("\n")
            epw.write(" dis_win_max = {} \n".format(self.epw_inputpara.dis_win_max))
            epw.write(" dis_froz_min = {}\n".format(self.epw_inputpara.dis_froz_min))
            epw.write(" dis_froz_max = {}\n".format(self.epw_inputpara.dis_froz_max))
            epw.write("\n")  
            
            epw.write(" fsthick  = {}             ! Fermi window [eV] : consider only states within Fermi energy +- fsthick\n".format(self.epw_inputpara.fsthick))
            epw.write(" degaussw      = {}        ! smearing in energy-conserving delta functions in [eV]\n".format(self.epw_inputpara.degaussw))
            epw.write(" degaussq      = {}        ! smearing for sum over q in the e-ph coupling in [meV]\n".format(self.epw_inputpara.degaussq))
            epw.write(" delta_qsmear  = {}        ! Change in the energy for each additional smearing in the a2f in [meV].\n".format(self.epw_inputpara.delta_qsmear))
            epw.write(" nqsmear       = {}        ! Number of different smearings used to calculate the a2f.\n".format(self.epw_inputpara.nqsmear))  
            epw.write(" phonselfen    = .true.    ! calculate the phonon self-energy\n")
            epw.write(" delta_approx  = .true.    ! apply the double delta approximation ! for phonon selfenergy.\n")
            epw.write(" nest_fn       = .true.    !calculate the nesting function.\n")
            epw.write(" selecqread    = .false.   ! If .true. then restart from the selecq.fmt file.\n")
            epw.write("                           ! selecq.fmt only needed if selecqread=.true. otherwise it will be re-created\n")
            epw.write("\n")

            epw.write(" iverbosity  = 2     ! 2 = verbose output for the superconducting part only.\n")
            epw.write("\n")
            epw.write(" nk1 = {}\n nk2 = {}\n nk3 = {}\n".format(self.epw_inputpara.nk[0], self.epw_inputpara.nk[1], self.epw_inputpara.nk[2]))
            epw.write("\n")
            epw.write(" nq1 = {}\n nq2 = {}\n nq3 = {}\n".format(self.epw_inputpara.nq[0], self.epw_inputpara.nq[1], self.epw_inputpara.nq[2]))
            epw.write("\n")

            epw.write(" nkf1 = {}\n nkf2 = {}\n nkf3 = {}\n".format(self.epw_inputpara.nkf[0], self.epw_inputpara.nkf[1], self.epw_inputpara.nkf[2]))
            epw.write("\n")
            epw.write(" filqf        = '{}'     \n".format(self.epw_inputpara.filqf))
            epw.write("! If any of elecselfen, phonselfen, specfun_el, or specfun_ph is true, mp_mesh_k must be false. The default value is false.\n")
            epw.write("/                           \n")                                                                

        return inputfilename
