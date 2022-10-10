import os
import logging
from argparse import ArgumentParser
from pathlib import Path
from pprint import pprint

from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.ase import AseAtomsAdaptor

from config import config
from specify_wyckoffs import specify_wyckoffs
from split_wyckoffs import split_wyckoffs
from substitution import substitution
from pso import pso
from psolib.utils.sort_atoms import sort_atoms

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
            while len(spe_wps.structs) < spe_wps.popsize:
                stru = spe_wps._gen_randomly()
                if hasattr(stru, "is_ordered"): # 判断结构是否是分数占据的无序结构
                    # 这里注意一定不要产生原胞结构，尽量保持晶胞结构
                    bstru = SpacegroupAnalyzer(stru).get_conventional_standard_structure()
                    spe_wps.structs.append(bstru)
                    logger.info(f"new you have successfully create No.{len(spe_wps.structs)} structures !")
            # write all the structures to the `work_path` by the format `.vasp` 
            if len(spe_wps.structs) == spe_wps.popsize:
                for i, struct in enumerate(spe_wps.structs):
                    logger.info(f"try successfully write POSCAR_{i+1} !")
                    filepath = os.path.join(spe_wps.work_path, "POSCAR_" + str(i+1))
                    _struct_ase = AseAtomsAdaptor.get_atoms(struct)
                    struct_ase = sort_atoms(_struct_ase, spe_wps.nameofatoms)
                    struct_ase.write(filepath, format='vasp')

        if self.config_d["mode"] == "splitwps":
            # init the parameter for generating the structure
            spl_wps = split_wyckoffs.init_from_config(self.config_d)
            # create structures and store them in the List `spe_wps.structs`
            while len(spl_wps.structs) < spl_wps.popsize:
                stru = spl_wps._gen_randomly()
                if hasattr(stru, "is_ordered"): # 判断结构是否是分数占据的无序结构
                    # 这里注意一定不要产生原胞结构，尽量保持晶胞结构
                    bstru = SpacegroupAnalyzer(stru).get_conventional_standard_structure()
                    spl_wps.structs.append(bstru)
                    logger.info(f"new you have successfully create No.{len(spl_wps.structs)} structures !")
            # write all the structures to the `work_path` by the format `.vasp` 
            if len(spl_wps.structs) == spl_wps.popsize:
                for i, struct in enumerate(spl_wps.structs):
                    logger.info(f"try successfully write POSCAR_{i+1} !")
                    filepath = os.path.join(spl_wps.work_path, "POSCAR_" + str(i+1))
                    _struct_ase = AseAtomsAdaptor.get_atoms(struct)
                    struct_ase = sort_atoms(_struct_ase, spl_wps.nameofatoms)
                    struct_ase.write(filepath, format='vasp')

        if self.config_d["mode"] == "substitution":
            substitution.init_from_config(self.config_d)

        if self.config_d["mode"] == "pso":
            if "specifywps" in self.sections:
                spe_wps = specify_wyckoffs.init_from_config(self.config_d)
                pso.init_from_config(self.config_d, specifywps=spe_wps)
            elif "splitwps" in self.sections:
                spl_wps = split_wyckoffs.init_from_config(self.config_d)
                pso.init_from_config(self.config_d, splitwps=spl_wps)