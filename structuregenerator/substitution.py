import time
from pathlib import Path
import logging

from itertools import product

import pandas as pd
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.vasp import Poscar


logger = logging.getLogger(__name__)

class substitution:

    def __init__(
        self,
        work_path,
        prototype_path,
        replacement,
        ):

        self.work_path = Path(work_path)
        self.prototype_path = prototype_path
        self.prototype_struct = Structure.from_file(self.prototype_path)
        # 获得原胞结构
        self.prototype_pstruct= SpacegroupAnalyzer(self.prototype_struct).get_primitive_standard_structure()
        self.replacement = replacement
    
        self.ternary_substitute()

    @classmethod
    def init_from_config(cls, config_d: dict):
        self = cls(
            work_path=config_d["work_path"],
            prototype_path=config_d["prototype_path"],
            replacement=config_d["replacement"],
        )
        return self
    
    def ternary_substitute(self):

        elements_list1 = self.replacement[0]
        elements_list2 = self.replacement[1]
        logger.info(f"The program will create {len(elements_list1)*len(elements_list2)} structutrs")
        time.sleep(3)
        num = 0
        serialnumber_formular = []
        for son_ele1 in elements_list1:
            parent_ele1 = elements_list1[0]
            for son_ele2 in elements_list2:
                num_formula = {}
                parent_ele2 = elements_list2[0]
                pstruct_copy = self.prototype_pstruct.copy()
                pstruct_copy.replace_species(
                    {
                        parent_ele1:son_ele1,
                        parent_ele2:son_ele2,
                        }
                    )
                num += 1
                logger.info(f"The program create No.{num}. {parent_ele1} -> {son_ele1} {parent_ele2} -> {son_ele2}")
                num_formula["serial_number"] = num
                num_formula["formula"] = "".join(pstruct_copy.composition.formula.split())
                serialnumber_formular.append(num_formula)
                #output_name = self.work_path.joinpath("POSCAR_"+str(num))
                output_name = self.work_path.joinpath(num_formula["formula"]+"_"+str(num)+".vasp")
                Poscar(pstruct_copy).write_file(output_name)
        
        struct_pd = pd.DataFrame(serialnumber_formular)
        struct_pd.to_csv(self.work_path.joinpath("structure.csv"), index=False)

