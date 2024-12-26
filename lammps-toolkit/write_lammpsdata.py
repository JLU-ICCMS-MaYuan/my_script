#!/usr/bin/env python3
import sys
import argparse

from ase.io import read
from ase.io.lammpsdata import write_lammps_data

parser = argparse.ArgumentParser(description='Generate supercells with optional diagonal or non-diagonal supercell matrices.')
# 输入文件
parser.add_argument('-i', '--input', type=str, default='POSCAR',
                        help='Input structure file (default: POSCAR)')

# 元素排序顺序
parser.add_argument('-o', '--order', nargs='+', type=str,
                    help='Custom element order for sorting, e.g., O C H for ordering oxygen first, followed by carbon and hydrogen, You have to make sure the order you specified is consistant with the order of input files')

args = parser.parse_args()

if "POSCAR" in args.input or ".vasp" in args.input:
    new_frame = read(args.input)
    write_lammps_data('data.lammps', new_frame, specorder=args.order)
elif ".lammpstrj" in args.input:
    # type_to_symbol = {i: j for i, j in enumerate(args.symbols)}  # i: [0,1,2], j: [Ce, Sr, H]
    # symbol_to_type = {v: k for k, v in type_to_symbol.items()}   # 用于输出cfg时的映射 v: [Ce, Sr, H], k: [0,1,2]
    atoms = read(open(args.input), format='lammps-dump-text', index=':')
    print(f"读取到 {len(atoms)} 帧轨迹")
    print("get last frame of trajectory")
    new_frame = atoms[-1]
    write_lammps_data('data.lammps', new_frame, specorder=args.order)