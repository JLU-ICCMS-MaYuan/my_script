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
            inputfilename = self.write_epw_phono_in(self.epw_inputpara.work_path)
            return inputfilename
        if mode == "epw_elph":
            inputfilename = self.write_epw_elph_in(self.epw_inputpara.work_path)
            return inputfilename
        if mode == "epw_aniso_sc":
            inputfilename = self.write_epw_aniso_sc_in(self.epw_inputpara.work_path)
            return inputfilename

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

            epw.write(" etf_mem     ={}\n".format(self.epw_inputpara.etf_mem))

            epw.write(" use_ws      = .false.\n")
            epw.write(" wannierize  = .true.\n")
            epw.write(" nbndsub     = {}\n".format(self.epw_inputpara.nbndsub))
            epw.write(" bands_skipped = 'exclude_bands = {}'\n".format(self.epw_inputpara.exclude_bands))
            epw.write(" num_iter    = {}\n".format(self.epw_inputpara.num_iter))
            epw.write(" dis_froz_min = {}\n".format(self.epw_inputpara.dis_froz_min))
            epw.write(" dis_froz_max = {}\n".format(self.epw_inputpara.dis_froz_max))
            for idx, pj in enumerate(self.epw_inputpara.proj):
                epw.write(" proj({})     = {}\n".format(idx+1, pj))
            epw.write(" wannier_plot= .true.\n")
            epw.write(" wdata(1)    = 'bands_plot = .true.'\n")
            epw.write(" wdata(2)    = 'begin kpoint_path'\n")
            for idx, path_name_coord in enumerate(self.epw_inputpara.path_name_coords_for_EPW):
                epw.write(f" wdata({idx+3})    = '{path_name_coord}'\n")
            epw.write(f" wdata({idx+4})    = 'end kpoint_path'\n")
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
            epw.write(" dvscf_dir   ='{}'".format(self.epw_inputpara.dvscf_dir))

            epw.write(" etf_mem     ={}\n".format(self.epw_inputpara.etf_mem))
            epw.write(" elph        = .true.\n")
            epw.write(" epbwrite    = .true.\n")
            epw.write(" epbread     = .false.\n")
            epw.write(" epwwrite    = .true.\n")
            epw.write(" epwread     = .false.\n")

            epw.write(" use_ws      = .false.\n")
            epw.write(" wannierize  = .true.\n")
            epw.write(" nbndsub     = {}\n".format(self.epw_inputpara.nbndsub))
            epw.write(" bands_skipped = 'exclude_bands = {}'\n".format(self.epw_inputpara.bands_skipped))
            epw.write(" num_iter    = {}\n".format(self.epw_inputpara.num_iter))
            epw.write(" dis_froz_min = {}\n".format(self.epw_inputpara.dis_froz_min))
            epw.write(" dis_froz_max = {}\n".format(self.epw_inputpara.dis_froz_max))
            for idx, pj in enumerate(self.epw_inputpara.proj):
                epw.write(" proj({})     = {}\n".format(idx, pj))
            epw.write(" wannier_plot= .true.\n")
            epw.write(" wdata(1)    = 'bands_plot = .true.'\n")
            epw.write(" wdata(2)    = 'begin kpoint_path'\n")
            for idx, path_name_coord in enumerate(self.epw_inputpara.path_name_coords_for_EPW):
                epw.write(f" wdata({idx+3}) = '{path_name_coord}'\n")
            epw.write(f" wdata({idx+4})     = 'end kpoint_path'\n")
            epw.write("/                           \n")      
        return inputfilename

    def write_epw_elph_in():
        pass

    def write_epw_aniso_sc_in():
        pass
