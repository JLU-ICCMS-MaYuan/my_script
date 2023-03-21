#!/usr/bin/env python3
import sys
from pathlib import Path

from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

print("使用格式:  symbolspg.py 目录路径")
print("例如:  symbolspg.py ./ ")
path = Path(input("请输入待计算的目录路径"))
for file in path.glob("*.vasp"):
    file_name = file.name
    struct = Structure.from_file(file_name)
    name = struct.composition.get_integer_formula_and_factor()[0]
    num = SpacegroupAnalyzer(struct).get_space_group_number()
    spg_symbol = SpacegroupAnalyzer(struct).get_space_group_symbol()
    print("{:<30}|  {:<5}|  {:<8}|  {:<10}".format(file_name, num, spg_symbol, name))