#!/public/home/mayuan/miniconda3/envs/cage/bin/python3
from ase.io.wien2k import write_struct
from ase.io.vasp import read_vasp
atoms = read_vasp("POSCAR")
write_struct("Ca.struct",atoms2=atoms)