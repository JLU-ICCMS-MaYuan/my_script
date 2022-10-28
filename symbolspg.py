#!/usr/bin/env python
import sys
from pathlib import Path

from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

path = Path(sys.argv[1])
for file in path.glob("U*"):
    file_name = file.name
    struct = Structure.from_file(file_name)
    name = struct.composition.get_integer_formula_and_factor()[0]
    num = SpacegroupAnalyzer(struct).get_space_group_number()
    spg_symbol = SpacegroupAnalyzer(struct).get_space_group_symbol()
    print("{:<20}|{:<5}|{:<8}|{:<10}".format(file_name, num, spg_symbol, name))