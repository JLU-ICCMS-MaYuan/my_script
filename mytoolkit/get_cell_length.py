#!/bin/bash

from pathlib import Path

from pymatgen.core.structure import Structure

for posfile in Path.cwd().rglob(f"POSCAR-*"):
    s = Structure.from_file("posfile")
    print("{} {}".format(s.formula.replace(' ', '') ,s.lattice.abs))
    