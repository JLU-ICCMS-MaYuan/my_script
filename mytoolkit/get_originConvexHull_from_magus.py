#!/usr/bin/env python3
import os
import sys

import pandas as pd

from pymatgen.analysis.phase_diagram import PDEntry, PhaseDiagram, PDPlotter
from pymatgen.core.composition import Composition
from pymatgen.core.periodic_table import Element

print("You can use it by:")
print("    python convert.py Ce Sc H 4 1 1")
print("    4 1 1 is hydrogen atoms number for CeH4 ScH H")
x, y, z = sys.argv[1:4]
a, b, c = sys.argv[4:]
a, b, c = float(a), float(b), float(c)

df = pd.read_table("extended_convex_hull", sep='\s+')

print("\n\n------------Atom as a endpoint for stable phase-----------")
with open("stable_by_atoms.dat", 'w') as ff:
    ff.write("{:<10} {:<4} {:<4} {:<4} {:<10}\n".format("formula", x, y, z, "enthalpy(eV/atom)"))
    print("{:<10} {:<4} {:<4} {:<4} {:<10}".format("formula", x, y, z, "enthalpy(eV/atom)"))
    for index, row in df.iterrows():
        formula = row['formula']
        comp = Composition(formula)
        ehull = row['ehull']
        x_atoms = comp[x]
        y_atoms = comp[y]
        z_atoms = comp[z]
        if ehull == 0 and len(comp)==3:
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

            ff.write("{:<10} {:<4} {:<4} {:<4}   {:<10.8f}\n".format(comp.formula.replace(" ", ""), x_atoms, y_atoms, z_atoms, ehull))
            print("{:<10} {:<4} {:<4} {:<4}   {:<10.8f}".format(comp.formula.replace(" ", ""), x_atoms, y_atoms, z_atoms, ehull))

print("\n\n------------Atom as a endpoint for unstable phase-----------")
with open("unstable_by_atoms.dat", 'w') as ff:
    ff.write("{:<10} {:<4} {:<4} {:<4} {:<10}\n".format("formula", x, y, z, "enthalpy(eV/atom)"))
    print("{:<10} {:<4} {:<4} {:<4} {:<10}".format("formula", x, y, z, "enthalpy(eV/atom)"))
    for index, row in df.iterrows():
        formula = row['formula']
        comp = Composition(formula)
        ehull = row['ehull']
        x_atoms = comp[x]
        y_atoms = comp[y]
        z_atoms = comp[z]
        if ehull > 0 and len(comp)==3:
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
                
            ff.write("{:<10} {:<4} {:<4} {:<4}   {:<10.8f}\n".format(comp.formula.replace(" ", ""), x_atoms, y_atoms, z_atoms, ehull))
            print("{:<10} {:<4} {:<4} {:<4}   {:<10.8f}".format(comp.formula.replace(" ", ""), x_atoms, y_atoms, z_atoms, ehull))




print("\n\n------------Block as a endpoint for stable phase-----------")
with open("stable_by_block.dat", 'w') as ff:
    ff.write("{:<10} {:<4} {:<4} {:<4} {:<10}\n".format("formula", x, y, z, "enthalpy(eV/atom)"))
    print("{:<10} {:<4} {:<4} {:<4} {:<10}".format("formula", x, y, z, "enthalpy(eV/atom)"))
    for index, row in df.iterrows():
        formula = row['formula']
        comp = Composition(formula)
        ehull = row['ehull']
        x_atoms = comp[x]
        y_atoms = comp[y]
        z_atoms = comp[z]
        if ehull == 0 and len(comp)==3:
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
                
            ff.write("{:<10} {:<4} {:<4} {:<4}   {:<10.8f}\n".format(comp.formula.replace(" ", ""), x_block, y_block, z_block, ehull))
            print("{:<10} {:<4} {:<4} {:<4}   {:<10.8f}".format(comp.formula.replace(" ", ""), x_block, y_block, z_block, ehull))

print("\n\n------------Block as a endpoint for unstable phase-----------")
with open("unstable_by_block.dat", 'w') as ff:
    ff.write("{:<10} {:<4} {:<4} {:<4} {:<10}\n".format("formula", x, y, z, "enthalpy(eV/atom)"))
    print("{:<10} {:<4} {:<4} {:<4} {:<10}".format("formula", x, y, z, "enthalpy(eV/atom)"))
    for index, row in df.iterrows():
        formula = row['formula']
        comp = Composition(formula)
        ehull = row['ehull']
        x_atoms = comp[x]
        y_atoms = comp[y]
        z_atoms = comp[z]
        if ehull > 0 and len(comp)==3:
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
                
            ff.write("{:<10} {:<4} {:<4} {:<4}   {:<10.8f}\n".format(comp.formula.replace(" ", ""), x_block, y_block, z_block, ehull))
            print("{:<10} {:<4} {:<4} {:<4}   {:<10.8f}".format(comp.formula.replace(" ", ""), x_block, y_block, z_block, ehull))



