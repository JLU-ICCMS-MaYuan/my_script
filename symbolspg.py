#!/usr/bin/env python3
import sys
from pathlib import Path

from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

print("使用格式:  symbolspg.py 目录路径 对称性的精度")
print("例如:  symbolspg.py ./ 0.01")

try:
    path = Path(sys.argv[1])
    symtol = eval(sys.argv[2])
except:
    path = Path(input("请输入待计算的目录路径\n"))
    symtol = 0.01

print("对称性的精度{}".format(str(symtol)))
for file in path.glob("POSCAR_*"):
    file_name = file.name
    struct = Structure.from_file(file_name)
    name = struct.composition.get_integer_formula_and_factor()[0]
    num = SpacegroupAnalyzer(struct, symprec=symtol).get_space_group_number()
    spg_symbol = SpacegroupAnalyzer(struct).get_space_group_symbol()
    print("{:<30}   {:<5}   {:<8}   {:<10}".format(file_name, num, spg_symbol, name))
