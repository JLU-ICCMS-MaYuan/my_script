#!/usr/bin/env python
import os
import sys
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

from ase.io import read, write

i = 0
for root, dirs, files in os.walk(os.getcwd()):
    for file in files:
        # print(file)
        if "cif" in file:
            atoms = read(os.path.join(root, file))
            numatom = atoms.get_number_of_atoms()
            formula = atoms.get_chemical_formula()
            print("{}  composition: {} numatom: {}".format(root, formula, numatom))
            # if numatom <= 25:
            i = i+1
            write("POSCAR_"+(str(i)), atoms, format="vasp")
