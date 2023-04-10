import os
import logging
from argparse import ArgumentParser

from pymatgen.io.ase import AseAtomsAdaptor

from structuregenerator.config import config
from structuregenerator.specify_wyckoffs import specify_wyckoffs 
from structuregenerator.split_wyckoffs import split_wyckoffs
from structuregenerator.split_wyckoffs_multprocessing import split_wyckoffs_multprocessing
from structuregenerator.substitution import substitution
from structuregenerator.pso import pso
from structuregenerator.psolib.utils.sort_atoms import sort_atoms

logger = logging.getLogger(__name__)

class generator_methods:

    def __init__(self, args: ArgumentParser) -> None:
        
        con = config(args)
        self.config_d = con.config_d
        self.sections = con.sections

        if self.config_d["mode"] == "specifywps":
            # init the parameter for generating the structure
            spe_wps = specify_wyckoffs.init_from_config(self.config_d)
            # create structures and store them in the List `spe_wps.structs`
            ################################# CASE 1 ##################################
            # CASE 1: create all the structures, then write them down !!!
            while len(spe_wps.structs) < spe_wps.popsize:
                stru = spe_wps._gen_randomly()
                if hasattr(stru, "is_ordered"): # 判断结构是否是分数占据的无序结构
                    _struct_ase = AseAtomsAdaptor.get_atoms(stru)
                    struct_ase = sort_atoms(_struct_ase, spe_wps.nameofatoms)
                    spe_wps.structs.append(struct_ase)
                    logger.info(f"new you have successfully create No.{len(spl_wps.structs)+1}-{str(_struct_ase.symbols)} !")
            # write all the structures to the `work_path` by the format `.vasp` 
            if len(spe_wps.structs) == spe_wps.popsize:
                for i, struct in enumerate(spe_wps.structs):
                    logger.info(f"try successfully write POSCAR_{i+1} !")
                    filepath = os.path.join(spe_wps.work_path, "POSCAR_" + str(i+1))
                    struct.write(filepath, format='vasp')
            ###########################################################################

        if self.config_d["mode"] == "splitwps":
            # init the parameter for generating the structure
            if self.config_d.get("max_workers", None):
                spl_wps = split_wyckoffs_multprocessing.init_from_config(self.config_d)
            else:
                spl_wps = split_wyckoffs.init_from_config(self.config_d)
            # create structures and store them in the List `spe_wps.structs`
            ################################# CASE 2 ##################################
            # CASE 2: create a structure, then write it !!!
            _struct_amounts = 0; clathrate = 0; ramdom_structure =0
            while _struct_amounts < spl_wps.popsize:
                stru, stru_type = spl_wps._gen_randomly()
                if hasattr(stru, "is_ordered"): # 判断结构是否是分数占据的无序结构
                    _struct_amounts += 1
                    if stru_type == "clathrate":
                        clathrate += 1
                    elif stru_type == "ramdom structure":
                        ramdom_structure += 1
                    _struct_ase = AseAtomsAdaptor.get_atoms(stru)
                    _struct_ase = sort_atoms(_struct_ase, spl_wps.nameofatoms)
                    _filepath = os.path.join(spl_wps.work_path, "POSCAR_" + str(_struct_amounts))
                    _struct_ase.write(_filepath, format='vasp')
                    logger.info(f"new you have successfully created No.{_struct_amounts}-{str(_struct_ase.symbols)}, and its type is {stru_type}!")
            ###########################################################################

        if self.config_d["mode"] == "substitution":
            substitution.init_from_config(self.config_d)

        if self.config_d["mode"] == "pso":
            if "specifywps" in self.sections:
                spe_wps = specify_wyckoffs.init_from_config(self.config_d)
                pso.init_from_config(self.config_d, specifywps=spe_wps)
            elif "splitwps" in self.sections:
                spl_wps = split_wyckoffs.init_from_config(self.config_d)
                pso.init_from_config(self.config_d, splitwps=spl_wps)