#!/usr/bin/env python3
import os
import sys
import shutil
from pathlib import Path
from ase.io import read

begin_id = int(sys.argv[1])
try:
    way_of_naming = sys.argv[2]
except:
    way_of_naming = 'retain_old_name'
    print("The default way of naming is 'retain_old_name'")

vaspfiles = list(Path.cwd().glob("*.vasp"))
vaspfiles = sorted(vaspfiles)

indexed_dirs = Path.cwd().joinpath("stlibs")
if not indexed_dirs.exists():
    os.mkdir(indexed_dirs) 

for id, old_filepath in enumerate(vaspfiles):
    s = read(old_filepath)
    current_id = id + begin_id
    old_file_name = old_filepath.name
    if way_of_naming == 'retain_old_name':
        new_file_name = str(current_id)+'.'+old_file_name.split('.')[0]+'.vasp'
    elif way_of_naming == 'only_formula':
        new_file_name = str(current_id)+'.'+s.get_chemical_formula()+'.vasp'
    else:
        new_file_name = str(current_id)+'.'+old_file_name.split('.')[0]+'.vasp'

    new_filepath  = indexed_dirs.joinpath(new_file_name)
    print(old_file_name, new_file_name)
    shutil.copy(old_filepath, new_filepath)
