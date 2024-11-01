#!/usr/bin/env python3
import sys
import argparse

from ase.io import read
from ase.build.tools import sort
from ase.io.lammpsdata import write_lammps_data

parser = argparse.ArgumentParser(description='Generate lammps inputed structure data')
# 输入文件
parser.add_argument('-i', '--input', type=str, default='POSCAR',
                        help='Input structure file (default: POSCAR)')

# 元素排序顺序
parser.add_argument('-o', '--order', nargs='+', type=str,
                    help='Custom element order for sorting, e.g., O C H for ordering oxygen first, followed by carbon and hydrogen, You have to make sure the order you specified is consistant with the order of input files')
args = parser.parse_args()

try:
    filename = sys.argv[1]
except:
    filename = 'POSCAR'

atom = read(args.input)
write_lammps_data('data.lammps', atom, specorder=args.order)
