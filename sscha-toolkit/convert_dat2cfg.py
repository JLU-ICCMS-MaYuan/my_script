#!/usr/bin/env python3
from pathlib import Path
import argparse
import numpy as np
import os
import sys

Ry2eV = 13.605698066
bohr2angstrom = 0.52917721092

def convert_dat2cfg(symbols, population):
    
    type_to_symbol = {i: j for i, j in enumerate(symbols)}  # i: [0,1,2], j: [Ce, Sr, H]
    symbol_to_type = {v: k for k, v in type_to_symbol.items()}   # 用于输出cfg时的映射 v: [Ce, Sr, H], k: [0,1,2]

    cfg_file = f"train_{str(population)}.cfg"
    cfg_file = open(cfg_file, "w")
    FORCE_FILES = True
    energy_file_name = Path("data_ensemble_manual").joinpath(f"energies_supercell_population{population}.dat")
    with open(energy_file_name, "r") as f:
        energy_lines = [l.strip() for l in f.readlines()]

    n_random = len(energy_lines)

    for j in range(1,n_random+1):
        structure_file_name = Path("data_ensemble_manual").joinpath(f"scf_population{population}_{j}.dat")
        with open(structure_file_name, "r") as f:
            lines = [l.strip() for l in f.readlines()]

        try:
            force_file_name = Path("data_ensemble_manual").joinpath(f"forces_population{population}_{j}.dat")
            with open(force_file_name, "r") as f:
                force_lines = [l.strip() for l in f.readlines()]
        except: 
            FORCE_FILES = False

        try:
            pressure_file_name = Path("data_ensemble_manual").joinpath(f"pressures_population{population}_{j}.dat")
            with open(pressure_file_name, "r") as f:
                pressure_lines = [l.strip() for l in f.readlines()]
        except:
            pressure_lines = [
            "0.0000 0.0000 0.0000",
            "0.0000 0.0000 0.0000",
            "0.0000 0.0000 0.0000",
            ]

        SIZE = len(lines)-6
        cell = np.zeros((3, 3))
        atoms = np.zeros((SIZE,3))
        forces = np.zeros((SIZE,3))
        atm_type = [None] * SIZE
        pressures = np.zeros((3, 3))

        for i in range(0,3):
            cell[i, :] = [float(x) for x in lines[i+1].split()[-3:]]

        for i in range(0,SIZE):
            atoms[i, :] = [float(x) for x in lines[i+6].split()[-3:]]
            atm_type[i] =  lines[i+6].split()[0]
            if FORCE_FILES == True:
                forces[i,:] = [float(x) for x in force_lines[i].split()[-3:]]
                forces[i,:] = forces[i, :] * Ry2eV / bohr2angstrom

        Volume = np.dot(cell[0],np.cross(cell[1], cell[2]))

        for i in range(0,3):

            pressures[i, :] = [float(x) for x in pressure_lines[i].split()[-3:]]
            pressures[i, :] = pressures[i, :] * Ry2eV / bohr2angstrom / bohr2angstrom / bohr2angstrom* Volume

        cfg_file.write("BEGIN_CFG\n")
        cfg_file.write(" Size\n")
        cfg_file.write("{: >5}".format(SIZE) +"\n")
        cfg_file.write(" Supercell\n")
        for row in cell:
            cfg_file.write("    {: >13f} {: >13f} {: >13f}\n".format(*row))
        cfg_file.write(" AtomData:  id type       cartes_x      cartes_y      cartes_z           fx          fy          fz\n")
        for i in range(0,len(atoms)):
            if FORCE_FILES == True:
                cfg_file.write("    {: >10}".format(i+1) + "{: >5}".format(symbol_to_type[atm_type[i]]) + "  {: >13f} {: >13f} {: >13f}".format(*atoms[i,:]) + "  {: >11f} {: >11f} {: >11f}\n".format(*forces[i,:]))
            else:
                cfg_file.write("    {: >10}".format(i+1) + "{: >5}".format(symbol_to_type[atm_type[i]]) + "  {: >13f} {: >13f} {: >13f}".format(*atoms[i,:]) + "  {: >11f} {: >11f} {: >11f}\n".format(*[0,0,0]))
        cfg_file.write(" Energy\n")
        cfg_file.write("     {: >13f}".format(float(energy_lines[j-1])*Ry2eV) +"\n")
        cfg_file.write(" PlusStress:  xx          yy          zz          yz          xz          xy\n")
        cfg_file.write("    {: >12f}".format(pressures[0,0]) + "{: >12f}".format(pressures[1,1]) +  "{: >12f}".format(pressures[2,2]))
        cfg_file.write("{: >12f}".format(pressures[1,2]) + "{: >12f}".format(pressures[0,2]) +  "{: >12f}\n".format(pressures[0,1]))
        cfg_file.write(" Feature atom_type_table " + str(symbol_to_type) + "\n")
        cfg_file.write(" Feature conf_number " + str(j) + "\n")
        cfg_file.write(" Feature population " + str(population) + "\n")
        cfg_file.write("END_CFG\n\n")
    cfg_file.close()
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-s","--symbols", type=str, nargs="+", default=["Ce", "Sc", "H"])
    parser.add_argument("-p", "--population", nargs="+", type=int, default=None)
    args = parser.parse_args()
    for population in args.population:
        convert_dat2cfg(args.symbols, population)
