#!/usr/bin/env python3
import os
import sys
from pathlib import Path
from ase.io import read
import numpy as np
from ase.atoms import Atoms
import argparse

def find_files(directory, pw_name_type=None, numberofpw=None):
    """
    查找目录中匹配特定模式的文件
    
    参数:
        directory (str): 搜索目录
        filename (str): 文件名匹配模式（支持前缀检测）
        numberofpw (int): 若文件名以 ESP_ 或 espresso_run_ 开头，生成指定数量的文件路径
    
    返回:
        list[Path]: 匹配的文件路径列表
    """
    result = []

    # 如果文件名以 ESP_ 开头且指定了数量
    if pw_name_type.startswith('ESP') and pw_name_type.endswith('pwi'):
        for num in range(numberofpw):
            pwfile = Path(directory).joinpath(f'ESP_{num}.pwi')
            result.append(pwfile)
        return result
    # 如果文件名以 espresso_run_ 开头且指定了数量
    elif pw_name_type.startswith('espresso_run') and pw_name_type.endswith('pwi'):
        for num in range(1, numberofpw+1):
            pwfile = Path(directory).joinpath(f'espresso_run_{num}.pwi')
            result.append(pwfile)
        return result
    elif pw_name_type.startswith('ESP') and pw_name_type.endswith('pwo'):
        for num in range(numberofpw):
            pwfile = Path(directory).joinpath(f'ESP_{num}.pwo')
            result.append(pwfile)
        return result
    # 如果文件名以 espresso_run_ 开头且指定了数量
    elif pw_name_type.startswith('espresso_run') and pw_name_type.endswith('pwo'):
        for num in range(1, numberofpw+1):
            pwfile = Path(directory).joinpath(f'espresso_run_{num}.pwo')
            result.append(pwfile)
        return result
    # 如果文件名以其他模式开头
    else:
        sys.exit(1)
        print("Error: 未提供有效的文件名")

def dump_cfg(frames, filename, symbol_to_type, mode='w'):
    with open(filename, mode) as f:
        for atoms in frames:
            ret = ''
            ret += 'BEGIN_CFG\n'
            ret += 'Size\n{: >5}\n'.format(len(atoms))
            try:
                cell = atoms.get_cell()[:]
                ret += 'Supercell\n{: >13f} {: >13f} {: >13f}\n{: >13f} {: >13f} {: >13f}\n{: >13f} {: >13f} {: >13f}\n'.format(*cell[0], *cell[1], *cell[2])
            except:
                pass
            cartes = atoms.positions
            has_forces = False
            has_constraints = False
            try:
                atoms.info['forces'] = atoms.get_forces()
            except:
                pass
            fields = ['id', 'type', 'cartes_x', 'cartes_y', 'cartes_z']
            if 'forces' in atoms.info:
                fields.extend(['fx', 'fy', 'fz'])
                forces = atoms.info['forces']
                has_forces = True

            flags = np.array(['T']*len(atoms))
            if atoms.constraints:
                atoms.info.append('mvable')
                for constr in atoms.constraints:
                    flags[constr.index] = 'F'
                has_constraints = True

            ret += 'AtomData: ' + ' '.join(fields) + '\n'
            for i, atom in enumerate(atoms):
                atom_info = '{: >10} {: >5} {: >13f} {: >13f} {: >13f} '.format(i + 1, symbol_to_type[atom.symbol], *cartes[i])
                if has_forces:
                    atom_info += '{: >11f} {: >11f} {: >11f} '.format(*forces[i])
                if has_constraints:
                    atom_info += flags[i]
                ret += atom_info + '\n'
            try:
                atoms.info['energy'] = atoms.get_potential_energy()
            except:
                pass
            if 'energy' in atoms.info:
                ret += 'Energy\n{: >13f}\n'.format(atoms.info['energy'])
            try:
                atoms.info['stress'] = atoms.get_stress()
            except:
                pass
            if 'stress' in atoms.info:
                stress = atoms.info['stress'] * atoms.get_volume() * -1.
                ret += 'PlusStress: xx yy zz yz xz xy\n{: >12f} {: >12f} {: >12f} {: >12f} {: >12f} {: >12f}\n'.format(*stress)
            if 'identification' in atoms.info:
                ret += 'Feature identification {}\n'.format(atoms.info['identification'])
            ret += 'END_CFG\n'
            f.write(ret)

    



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Convert QE output to cfg file')
    parser.add_argument('-r', '--root_dir',  default='', help='root directory of QE output')
    parser.add_argument('-n', '--number_of_files', default=None, type=int, help='number of QE output files')
    parser.add_argument('-s', '--symbols',   default='', nargs='+', help='name of elements')
    parser.add_argument('-p', '--pw_name_type', default='', help='ESP.pwo, ESP.pwi, espresso_run.pwo, espresso_run.pwi')
    args = parser.parse_args()
    type_to_symbol = {i: j for i, j in enumerate(args.symbols)}  # i: [0,1,2], j: [Ce, Sr, H]
    symbol_to_type = {v: k for k, v in type_to_symbol.items()}   # 用于输出cfg时的映射 v: [Ce, Sr, H], k: [0,1,2]
    pwfiles = find_files(args.root_dir, args.pw_name_type, args.number_of_files)
    frames  = [read(pwfile) for pwfile in pwfiles]
    # print(frames)
    dump_cfg(frames, 'train.cfg', symbol_to_type, mode='w')