#!/usr/bin/env python3
import sys

from ase.io import read
from ase.io.lammpsdata import write_lammps_data

try:
    filename = sys.argv[1]
except:
    filename = 'POSCAR'

atom = read(filename)
write_lammps_data(atom)