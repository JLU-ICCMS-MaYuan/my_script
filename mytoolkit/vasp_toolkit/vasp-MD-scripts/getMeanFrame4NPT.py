import sys

from ase.io.vasp import read_vasp_xdatcar, write_vasp 
from ase import Atoms 
import numpy as np 

try:
    filename = sys.argv[1]
except:
    filename = "XDATCAR"
    
begin_steps = int(sys.argv[2])
end_steps   = int(sys.argv[3])

frames = read_vasp_xdatcar(filename, index=slice(begin_steps, end_steps))

position_list = []
cell_list = []
for frame in frames:
    position_list.append(frame.positions @ np.linalg.inv(frame.cell[:]))
    cell_list.append(frame.cell[:])

position_mean = np.mean(position_list, axis=0)
cell_mean = np.mean(cell_list, axis=0)
symbols = frame.symbols 

frame_mean = Atoms(symbols=symbols,
                   positions=position_mean @ cell_mean,
                   cell=cell_mean)

write_vasp("npt-mean.vasp", frame_mean)