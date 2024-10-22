#!/usr/bin/env python3
import sys

from ase.io.vasp import read_vasp_xdatcar, write_vasp 
from ase import Atoms 
import numpy as np

try:
    filename = sys.argv[1]
except:
    filename = "XDATCAR"
    
begin_steps = int(sys.argv[2]) # 包含左
end_steps   = int(sys.argv[3]) # 不包含后右

frames = read_vasp_xdatcar(filename, index=slice(begin_steps, end_steps))

position_list = []
for frame in frames:
    position_list.append(frame.positions)

position_mean = np.mean(position_list, axis=0)
cell = frame.cell[:]
symbols = frame.symbols 

frame_mean = Atoms(symbols=symbols,
                   positions=position_mean,
                   cell=cell)

write_vasp("nvt-mean.vasp", frame_mean)
