import os
import sys

from pymatgen.analysis.phase_diagram import PDEntry, PhaseDiagram, PDPlotter
from pymatgen.core.composition import Composition
from pymatgen.core.periodic_table import Element

os.system("cp extended_convex_hull extended_convex_hull.dat")
os.system("sed '1,6d' extended_convex_hull.dat")
x, y, z = sys.argv[1:]

with open("extended_convex_hull.dat", 'r') as f:
    lines = f.readlines()

with open("stable.dat", 'w') as ff:
    ff.write("{:<10} {:<4} {:<4} {:<4} {:<10}\n".format("formula", x, y, z, "enthalpy(eV/atom)"))
    print("{:<10} {:<4} {:<4} {:<4} {:<10}".format("formula", x, y, z, "enthalpy(eV/atom)"))
    for line in lines:
        #print(x+line.split()[2] + y+line.split()[3] + z+line.split()[4])
        formula = x+line.split()[2] + y+line.split()[3] + z+line.split()[4]
        fitness = float(line.split()[8])
        comp = Composition(formula)
        if fitness == 0 and len(comp.element_composition)==3:
            try:
                x_atoms = comp.to_data_dict['unit_cell_composition'][x]
            except:
                x_atoms = 0

            try:
                y_atoms = comp.to_data_dict['unit_cell_composition'][y]
            except:
                y_atoms = 0

            try:
                z_atoms = comp.to_data_dict['unit_cell_composition'][z]
            except:
                z_atoms = 0

            ff.write("{:<10} {:<4} {:<4} {:<4}   {:<10.8f}\n".format(comp.formula.replace(" ", ""), x_atoms, y_atoms, z_atoms, fitness))
            print("{:<10} {:<4} {:<4} {:<4}   {:<10.8f}".format(comp.formula.replace(" ", ""), x_atoms, y_atoms, z_atoms, fitness))

with open("unstable.dat", 'w') as ff:
    ff.write("{:<10} {:<4} {:<4} {:<4} {:<10}\n".format("formula", x, y, z, "enthalpy(eV/atom)"))
    print("{:<10} {:<4} {:<4} {:<4} {:<10}".format("formula", x, y, z, "enthalpy(eV/atom)"))
    for line in lines:
        #print(x+line.split()[2] + y+line.split()[3] + z+line.split()[4])
        formula = x+line.split()[2] + y+line.split()[3] + z+line.split()[4]
        fitness = float(line.split()[8])
        comp = Composition(formula)
        if fitness > 0 and len(comp.element_composition)==3:
            try:
                x_atoms = comp.to_data_dict['unit_cell_composition'][x]
            except:
                x_atoms = 0

            try:
                y_atoms = comp.to_data_dict['unit_cell_composition'][y]
            except:
                y_atoms = 0

            try:
                z_atoms = comp.to_data_dict['unit_cell_composition'][z]
            except:
                z_atoms = 0
                
            ff.write("{:<10} {:<4} {:<4} {:<4}   {:<10.8f}\n".format(comp.formula.replace(" ", ""), x_atoms, y_atoms, z_atoms, fitness))
            print("{:<10} {:<4} {:<4} {:<4}   {:<10.8f}".format(comp.formula.replace(" ", ""), x_atoms, y_atoms, z_atoms, fitness))


with open("stable_by_block.dat", 'w') as ff:
    ff.write("{:<10} {:<4} {:<4} {:<4} {:<10}\n".format("formula", x, y, z, "enthalpy(eV/atom)"))
    print("{:<10} {:<4} {:<4} {:<4} {:<10}".format("formula", x, y, z, "enthalpy(eV/atom)"))
    for line in lines:
        #print(x+line.split()[2] + y+line.split()[3] + z+line.split()[4])
        formula = x+line.split()[2] + y+line.split()[3] + z+line.split()[4]
        fitness = float(line.split()[8])
        comp = Composition(formula)
        if fitness == 0 and len(comp.element_composition)==3:
            try:
                x_atoms = comp.to_data_dict['unit_cell_composition'][x]
                x_block = x_atoms
            except:
                x_block = 0

            try:
                y_atoms = comp.to_data_dict['unit_cell_composition'][y]
                y_block = y_atoms
            except:
                y_block = 0

            try:
                z_atoms = comp.to_data_dict['unit_cell_composition'][z]
                z_block = z_atoms - 4*x_block - 2*y_block
                if z_block < 0:
                    break
            except:
                z_block = 0
                
            ff.write("{:<10} {:<4} {:<4} {:<4}   {:<10.8f}\n".format(comp.formula.replace(" ", ""), x_block, y_block, z_block, fitness))
            print("{:<10} {:<4} {:<4} {:<4}   {:<10.8f}".format(comp.formula.replace(" ", ""), x_block, y_block, z_block, fitness))

with open("unstable_by_block.dat", 'w') as ff:
    ff.write("{:<10} {:<4} {:<4} {:<4} {:<10}\n".format("formula", x, y, z, "enthalpy(eV/atom)"))
    print("{:<10} {:<4} {:<4} {:<4} {:<10}".format("formula", x, y, z, "enthalpy(eV/atom)"))
    for line in lines:
        #print(x+line.split()[2] + y+line.split()[3] + z+line.split()[4])
        formula = x+line.split()[2] + y+line.split()[3] + z+line.split()[4]
        fitness = float(line.split()[8])
        comp = Composition(formula)
        if fitness > 0 and len(comp.element_composition)==3:
            try:
                x_atoms = comp.to_data_dict['unit_cell_composition'][x]
                x_block = x_atoms
            except:
                x_block = 0
                
            try:
                y_atoms = comp.to_data_dict['unit_cell_composition'][y]
                y_block = y_atoms
            except:
                y_block = 0

            try:
                z_atoms = comp.to_data_dict['unit_cell_composition'][z]
                z_block = z_atoms - 4*x_block - 2*y_block
                if z_block < 0:
                    break
            except:
                z_block = 0
                
            ff.write("{:<10} {:<4} {:<4} {:<4}   {:<10.8f}\n".format(comp.formula.replace(" ", ""), x_block, y_block, z_block, fitness))
            print("{:<10} {:<4} {:<4} {:<4}   {:<10.8f}".format(comp.formula.replace(" ", ""), x_block, y_block, z_block, fitness))


