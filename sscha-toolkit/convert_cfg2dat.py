#!/usr/bin/env python3
import os
import numpy as np

import argparse

from collections import defaultdict
from ase import Atoms
from ase.units import GPa


def load_cfg(filename, type_to_symbol):
    frames = []
    with open(filename) as f:
        line = 'chongchongchong!'
        while line:
            line = f.readline()
            if 'BEGIN_CFG' in line:
                cell = np.zeros((3, 3))

            if 'Size' in line:
                line = f.readline()
                natoms = int(line.split()[0])
                positions = np.zeros((natoms, 3))
                forces = np.zeros((natoms, 3))
                energies = np.zeros(natoms)
                symbols = ['X'] * natoms
                constraints = ['T'] * natoms

            if 'Supercell' in line: 
                for i in range(3):
                    line = f.readline()
                    for j in range(3):
                        cell[i, j] = float(line.split()[j])

            if 'AtomData' in line:
                d = defaultdict(int)
                for (i, x) in enumerate(line.split()[1:]):
                    d[x] = i

                for _ in range(natoms):
                    line = f.readline()
                    fields = line.split()
                    i = int(fields[d['id']]) - 1
                    symbols[i] = type_to_symbol[int(fields[d['type']])]
                    positions[i] = [float(fields[d[attr]]) for attr in ['cartes_x', 'cartes_y' ,'cartes_z']]
                    forces[i] = [float(fields[d[attr]]) for attr in ['fx', 'fy' ,'fz']]
                    energies[i] = float(fields[d['site_en']])

                atoms = Atoms(symbols=symbols, cell=cell, positions=positions, pbc=True)
                if d['fx'] != 0:
                    atoms.info['forces'] = forces
                if d['site_en'] != 0:
                    atoms.info['energies'] = energies

            if 'Energy' in line and 'Weight' not in line:
                line = f.readline()
                atoms.info['energy'] = float(line.split()[0])

            if 'PlusStress' in line:
                line = f.readline()
                plus_stress = np.array(list(map(float, line.split())))
                atoms.info['PlusStress'] = plus_stress
                atoms.info['stress'] = -plus_stress / atoms.get_volume()
                atoms.info['pstress'] = atoms.info['stress'] / GPa

            if 'END_CFG' in line:
                frames.append(atoms)

            if 'identification' in line:
                atoms.info['identification'] = int(line.split()[2])

    return frames


def convert_cfg2dat(frames, ensembles_name, sscha_pop):

    angstrom2bohr = 1.88972612462577
    eV2Ry = 0.0734986443513116
    
    energies = np.zeros(len(frames))
    

    for i, atoms in enumerate(frames):
        energies[i] = atoms.info['energy']*eV2Ry
        force_filename = os.path.join(ensembles_name, f"forces_population{sscha_pop}_{i}.dat")
        np.savetxt(force_filename, atoms.info['forces']*(eV2Ry/angstrom2bohr))
        
        pressure_filename = os.path.join(ensembles_name, f"pressures_population{sscha_pop}_{i}.dat")
        pressure_eV_in_a_volume_Angstrom3 = atoms.info['PlusStress']
        pressure_Ry_per_volume_bohr3 = [[0,0,0],[0,0,0],[0,0,0]]
        pressure_Ry_per_volume_bohr3[0][0] = pressure_eV_in_a_volume_Angstrom3[0]  *  eV2Ry/(atoms.get_volume()*angstrom2bohr**3)
        pressure_Ry_per_volume_bohr3[1][1] = pressure_eV_in_a_volume_Angstrom3[1]  *  eV2Ry/(atoms.get_volume()*angstrom2bohr**3)     
        pressure_Ry_per_volume_bohr3[2][2] = pressure_eV_in_a_volume_Angstrom3[2]  *  eV2Ry/(atoms.get_volume()*angstrom2bohr**3)
        pressure_Ry_per_volume_bohr3[1][2] = pressure_eV_in_a_volume_Angstrom3[3]  *  eV2Ry/(atoms.get_volume()*angstrom2bohr**3)     
        pressure_Ry_per_volume_bohr3[0][2] = pressure_eV_in_a_volume_Angstrom3[4]  *  eV2Ry/(atoms.get_volume()*angstrom2bohr**3)
        pressure_Ry_per_volume_bohr3[0][1] = pressure_eV_in_a_volume_Angstrom3[5]  *  eV2Ry/(atoms.get_volume()*angstrom2bohr**3)     
        np.savetxt(pressure_filename, pressure_Ry_per_volume_bohr3)

    energy_file_name = os.path.join(ensembles_name, f"energies_supercell_population{sscha_pop}.dat")
    np.savetxt(energy_file_name, energies)
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Convert cfg to dat')
    parser.add_argument('-i', '--input-cfg',  default='', help='input cfg file')
    parser.add_argument('-s', '--symbols',    default='', nargs='+', help='name of elements')
    parser.add_argument('-p', '--population',  default='', type=int,  help='sscha population')
    parser.add_argument('-e', '--ensembles',  default='ensembles', help='ensembles name')
    args = parser.parse_args()
    type_to_symbol = {i: j for i, j in enumerate(args.symbols)}  # i: [0,1,2], j: [Ce, Sr, H]
    symbol_to_type = {v: k for k, v in type_to_symbol.items()}   # 用于输出cfg时的映射 v: [Ce, Sr, H], k: [0,1,2]

    frames = load_cfg(args.input_cfg, type_to_symbol)
    if not os.path.exists(args.ensembles):
        os.mkdir(args.ensembles)
    convert_cfg2dat(frames, args.ensembles, args.population)
