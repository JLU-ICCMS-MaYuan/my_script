#!/usr/bin/env python
import os
import sys
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

from pathlib import Path
from ase.io import read, write

i = 0
for root, dirs, files in os.walk(os.getcwd()):
    # print(root)
    if str(100) in root:
        for file in files:
            if "POSCAR" ==  file:
                atoms = read(os.path.join(root, file))
                numatom = atoms.get_global_number_of_atoms()
                formula = atoms.get_chemical_formula()
                # print("{}  composition: {} numatom: {}".format(root, formula, numatom))
                print("{},{},{}".format(Path(root).parent.name.split('-')[1], formula, numatom))
                # if numatom <= 25:
                i = i+1
            else:
                print(root)
            # write("POSCAR_"+(str(i)), atoms, format="vasp")


