#!/usr/bin/env python3
import sys

from pymatgen.analysis.phase_diagram import PDEntry, PhaseDiagram, PDPlotter
from pymatgen.core.composition import Composition
from pymatgen.core.periodic_table import Element

print("You can use it by:")
print("    python convert.py Ce Sr H 4 2 1")
x, y, z = sys.argv[1:4]
a, b, c = sys.argv[4:]
a, b, c = float(a), float(b), float(c)

with open("stable.csv", 'r') as f:
    lines = f.readlines()

print("\n\n------------Atom as a endpoint for stable phase-----------")
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

print("\n\n------------Block as a endpoint for stable phase-----------")
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
                z_block = z_atoms - a*x_block - b*y_block
                if z_block < 0:
                    continue
            except:
                z_block = 0
                
            ff.write("{:<8} {:<10} {:<4} {:<4} {:<4}   {:<10.8f}\n".format(idx, comp.formula.replace(" ", ""), x_block, y_block, z_block, fitness))
            print("{:<8} {:<10} {:<4} {:<4} {:<4}   {:<10.8f}".format(idx, comp.formula.replace(" ", ""), x_block, y_block, z_block, fitness))



with open("unstable.csv", 'r') as f:
    lines = f.readlines()

print("\n\n------------Atom as a endpoint for unstable phase-----------")
with open("unstable_by_atoms.dat", 'w') as ff:
    ff.write("{:<8} {:<10} {:<4} {:<4} {:<4} {:<10}\n".format("idx", "formula", x, y, z, "enthalpy(eV/atom)"))
    print("{:<8} {:<10} {:<4} {:<4} {:<4} {:<10}".format("idx", "formula", x, y, z, "enthalpy(eV/atom)"))
    for line in lines:
        idx     = line.split(',')[0]
        formula = line.split(',')[1]
        fitness = float(line.split(',')[-1])
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
                
            ff.write("{:<8} {:<10} {:<4} {:<4} {:<4}   {:<10.8f}\n".format(idx, comp.formula.replace(" ", ""), x_atoms, y_atoms, z_atoms, fitness))
            print("{:<8} {:<10} {:<4} {:<4} {:<4}   {:<10.8f}".format(idx, comp.formula.replace(" ", ""), x_atoms, y_atoms, z_atoms, fitness))

print("\n\n------------Block as a endpoint for unstable phase-----------")
with open("unstable_by_block.dat", 'w') as ff:
    ff.write("{:<8} {:<10} {:<4} {:<4} {:<4} {:<10}\n".format("idx", "formula", x, y, z, "enthalpy(eV/atom)"))
    print("{:<8} {:<10} {:<4} {:<4} {:<4} {:<10}".format("idx", "formula", x, y, z, "enthalpy(eV/atom)"))
    for line in lines:
        idx     = line.split(',')[0]
        formula = line.split(',')[1]
        fitness = float(line.split(',')[-1])
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
                z_block = z_atoms - a*x_block - b*y_block
                if z_block < 0:
                    continue
            except:
                z_block = 0
                
            ff.write("{:<8} {:<10} {:<4} {:<4} {:<4}   {:<10.8f}\n".format(idx, comp.formula.replace(" ", ""), x_block, y_block, z_block, fitness))
            print("{:<8} {:<10} {:<4} {:<4} {:<4}   {:<10.8f}".format(idx, comp.formula.replace(" ", ""), x_block, y_block, z_block, fitness))


