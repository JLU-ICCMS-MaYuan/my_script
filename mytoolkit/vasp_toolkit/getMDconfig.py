#!/usr/bin/env python3
import sys

from ase.io.vasp import read_vasp_xdatcar, write_vasp 

try:
    filename = sys.argv[1]
except:
    filename = "XDATCAR"

steps = int(sys.argv[2])

frames = read_vasp_xdatcar(filename, index=steps)

write_vasp("POSCAR_"+str(steps)+".vasp", frames)