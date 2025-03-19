#!/usr/bin/env python3

import argparse

import numpy as np
from ase.atoms import Atoms
from ase.units import GPa
from collections import defaultdict

import spglib

#TODO 
# different cell
# pbc
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
       
if __name__ == "__main__":
    
    parser = argparse.ArgumentParser("read_cfg.py -dft train.cfg  -mlp calculated.cfg   -s Ce Sr H -t 500 10 -o train1.cfg")
    parser.add_argument('-i', '--input', default='new_train.cfg', help='input cfg filename')
    parser.add_argument('-s', '--symbols', default='', nargs='+', help='name of elements')
    parser.add_argument('-p', '--prec', default=0.01, type=float, help="get tolerance of symmetry ")

    args = parser.parse_args()
    
    type_to_symbol = {i: j for i, j in enumerate(args.symbols)}  # i: [0,1,2], j: [Ce, Sr, H]
    symbol_to_type = {v: k for k, v in type_to_symbol.items()}   # 用于输出cfg时的映射 v: [Ce, Sr, H], k: [0,1,2]

    # 载入数据, 并判断dftcfg的结构数是否等于mlpcfg的结构数
    DFTframes = load_cfg(args.input, type_to_symbol)

    for idx, frame in enumerate(DFTframes):
        # 获取原子总数
        total_atoms = len(frame)
        # 获取空间群对称性
        lattice = frame.get_cell()
        positions = frame.get_scaled_positions()
        numbers = frame.get_atomic_numbers()
        cell = (lattice, positions, numbers)
        spacegroup = spglib.get_spacegroup(cell, args.prec)
        pstress = -sum(frame.info['pstress'][0:3])/3
        # 打印原子数和空间群对称性
        print("{:<10} {:<5} {:<15}  {:>6.3f}  {:>6.2f}".format(str(frame.symbols), total_atoms, spacegroup, float(frame.get_volume()), pstress))
