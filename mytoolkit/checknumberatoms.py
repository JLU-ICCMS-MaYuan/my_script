#!/usr/bin/env python
import os
import sys

from pymatgen.core.structure import Structure

for root, dirs, files in os.walk(os.getcwd()):
    if "POSCAR" in files:
        struct = Structure.from_file(os.path.join(root, "POSCAR"))
        comp = dict(struct.composition.get_el_amt_dict())
        numatom = sum(list(comp.values()))
        print("{}  composition: {} numatom: {}".format(root, comp, numatom))
