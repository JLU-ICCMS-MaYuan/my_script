#!/usr/bin/env python3
import sys
import os

from pathlib import Path

from ase.io.vasp import read_vasp_xdatcar, write_vasp 

dirnames = sys.argv[1:]

dirnames = ['2000K/1.potim_0.1_and_smass_2/', '2000K/2.potim_0.1_and_smass_0/', '2000K/3.potim_0.1_and_smass_0/', '2500/', '2500/2.potim_0.1_and_smass_0/', '3000K/1.potim_0.1_and_smass_2/', '3000K/2.potim_0.1_and_smass_0/', '3500/' '4000K/1.potim_0.1_and_smass_2/', '4000K/2.potim_0.1_and_smass_0/']

if not os.path.exists('slices.poscar'):
    os.mkdir('slices.poscar')

index = 0
for dirname in dirnames:
    print(dirname)
    xdatcar_path = Path(dirname).joinpath("XDATCAR")
    num_configs = os.popen("grep 'Direct configuration' {} | wc -l".format(xdatcar_path)).read().strip()
    print('Total steps: {}'.format(num_configs))
    xdatcar = read_vasp_xdatcar(xdatcar_path, slice(0, int(num_configs)))
    for step in range(0, num_configs, 5):
        poscar = xdatcar[step]
        index += 1
        write_vasp("slices.poscar/"+ "POSCAR_"+str(index)+".vasp", poscar)
        print()