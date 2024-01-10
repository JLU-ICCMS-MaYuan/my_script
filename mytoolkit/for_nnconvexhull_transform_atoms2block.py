import sys

from pymatgen.analysis.phase_diagram import PDEntry, PhaseDiagram, PDPlotter
from pymatgen.core.composition import Composition
from pymatgen.core.periodic_table import Element

x, y, z = sys.argv[1:]

with open("stable.csv", 'r') as f:
    lines = f.readlines()

with open("unstable.csv", 'r') as f:
    lines = f.readlines()

with open("stable_by_atoms.dat", 'w') as ff:
    ff.write("{:<8} {:<10} {:<4} {:<4} {:<4} {:<10}\n".format("idx", "formula", x, y, z, "enthalpy(eV/atom)"))
    print("{:<8} {:<10} {:<4} {:<4} {:<4} {:<10}".format("idx", "formula", x, y, z, "enthalpy(eV/atom)"))
    for line in lines:
        idx     = line.split(',')[0]
        formula = line.split(',')[1]
        fitness = float(line.split(',')[-1])
        comp = Composition(formula)
        if fitness == 0:
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
            ff.write("{:<8} {:<10} {:<4} {:<4} {:<4}   {:<10.8f}\n".format(idx, comp.formula.replace(" ", ""), x_atoms, y_atoms, z_atoms, fitness))
            print("{:<8} {:<10} {:<4} {:<4} {:<4}   {:<10.8f}".format(idx, comp.formula.replace(" ", ""), x_atoms, y_atoms, z_atoms, fitness))


with open("unstable_by_atoms.dat", 'w') as ff:
    ff.write("{:<8} {:<10} {:<4} {:<4} {:<4} {:<10}\n".format("idx", "formula", x, y, z, "enthalpy(eV/atom)"))
    print("{:<8} {:<10} {:<4} {:<4} {:<4} {:<10}".format("idx", "formula", x, y, z, "enthalpy(eV/atom)"))
    for line in lines:
        idx     = line.split(',')[0]
        formula = line.split(',')[1]
        fitness = float(line.split(',')[-1])
        comp = Composition(formula)
        if fitness > 0:
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
                
            ff.write("{:<8} {:<10} {:<4} {:<4} {:<4}   {:<10.8f}\n".format(idx, comp.formula.replace(" ", ""), x_atoms, y_atoms, z_atoms, fitness))
            print("{:<8} {:<10} {:<4} {:<4} {:<4}   {:<10.8f}".format(idx, comp.formula.replace(" ", ""), x_atoms, y_atoms, z_atoms, fitness))


with open("stable_by_block.dat", 'w') as ff:
    ff.write("{:<8} {:<10} {:<4} {:<4} {:<4} {:<10}\n".format("idx", "formula", x, y, z, "enthalpy(eV/atom)"))
    print("{:<8} {:<10} {:<4} {:<4} {:<4} {:<10}".format("idx", "formula", x, y, z, "enthalpy(eV/atom)"))
    for line in lines:
        idx     = line.split(',')[0]
        formula = line.split(',')[1]
        fitness = float(line.split(',')[-1])
        comp = Composition(formula)
        if fitness == 0:
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
            except:
                z_block = 0
                
            ff.write("{:<8} {:<10} {:<4} {:<4} {:<4}   {:<10.8f}\n".format(idx, comp.formula.replace(" ", ""), x_block, y_block, z_block, fitness))
            print("{:<8} {:<10} {:<4} {:<4} {:<4}   {:<10.8f}".format(idx, comp.formula.replace(" ", ""), x_block, y_block, z_block, fitness))


with open("unstable_by_block.dat", 'w') as ff:
    ff.write("{:<8} {:<10} {:<4} {:<4} {:<4} {:<10}\n".format("idx", "formula", x, y, z, "enthalpy(eV/atom)"))
    print("{:<8} {:<10} {:<4} {:<4} {:<4} {:<10}".format("idx", "formula", x, y, z, "enthalpy(eV/atom)"))
    for line in lines:
        idx     = line.split(',')[0]
        formula = line.split(',')[1]
        fitness = float(line.split(',')[-1])
        comp = Composition(formula)
        if fitness > 0:
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
            except:
                z_block = 0
                
            ff.write("{:<8} {:<10} {:<4} {:<4} {:<4}   {:<10.8f}\n".format(idx, comp.formula.replace(" ", ""), x_block, y_block, z_block, fitness))
            print("{:<8} {:<10} {:<4} {:<4} {:<4}   {:<10.8f}".format(idx, comp.formula.replace(" ", ""), x_block, y_block, z_block, fitness))


