#!/usr/bin/env python3
import sys
import argparse

from ase.io import read
from ase.build.tools import sort
from ase.io.lammpsdata import write_lammps_data
from ase.io.lammpsrun import read_lammps_dump_text

parser = argparse.ArgumentParser(description='Generate supercells with optional diagonal or non-diagonal supercell matrices.')
# 输入文件
parser.add_argument('-i', '--input', type=str, default='POSCAR',
                        help='Input structure file (default: POSCAR)')

# 元素排序顺序
parser.add_argument('-o', '--order', nargs='+', type=str,
                    help='Custom element order for sorting, e.g., O C H for ordering oxygen first, followed by carbon and hydrogen, You have to make sure the order you specified is consistant with the order of input files')
args = parser.parse_args()

if "POSCAR" in args.input or ".vasp" in args.input:
    atoms = read(args.input)
elif ".lammpstrj" in args.input:
    atoms = read_lammps_dump_text(open(args.input), index=-1)
    print("get last frame of trajectory")
    
write_lammps_data('data.lammps', atoms, specorder=args.order)
