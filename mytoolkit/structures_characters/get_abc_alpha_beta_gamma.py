#!/usr/bin/env python3

import sys
from pymatgen.core.structure import Structure

# 读取结构文件路径
filename = sys.argv[1]

# 读取结构
structure = Structure.from_file(filename)

# 获取晶格信息
lattice = structure.lattice

# 输出晶胞参数
print(f"a = {lattice.a:.6f} Å")
print(f"b = {lattice.b:.6f} Å")
print(f"c = {lattice.c:.6f} Å")
print(f"alpha = {lattice.alpha:.6f} °")
print(f"beta  = {lattice.beta:.6f} °")
print(f"gamma = {lattice.gamma:.6f} °")
print(f"volume = {lattice.volume}")