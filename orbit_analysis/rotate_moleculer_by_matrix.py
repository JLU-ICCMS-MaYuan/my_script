#!/usr/bin/env python3
import numpy as np

import sys

xyzfile = sys.argv[1]
# 旋转矩阵：请按需替换成你要的矩阵
# rotation_matrix = np.array([
#     [1., 0., 0.],
#     [0., 0.577228, -0.816583],
#     [0., 0.816583, 0.577228]
# ])
rotation_matrix = np.array([
    [0.70710, -0.70710, 0],
    [0.70710,  0.70710, 0],
    [0,        0,       1]
])
with open(xyzfile, "r") as f:
    coords_str = f.read()

# 解析坐标
atoms = []
for line in coords_str.strip().splitlines():
    parts = line.split()
    if len(parts) == 4:
        symbol = parts[0]
        pos = list(map(float, parts[1:]))
        atoms.append((symbol, np.array(pos)))

# 应用旋转
rotated_atoms = []
index_atoms = []
symbols_atoms = []
for i, (symbol, pos) in enumerate(atoms, start=1):
    new_pos = rotation_matrix @ pos
    index_atoms.append(f"A{i}")
    rotated_atoms.append(new_pos)
    symbols_atoms.append(symbol)

# 打印结果
for index, vec in zip(index_atoms, rotated_atoms):
    print(f"{index} = {{{vec[0]: .6f}, {vec[1]: .6f}, {vec[2]: .6f}}};")

for symbol, vec in zip(symbols_atoms, rotated_atoms):
    print(f"{symbol:<2}  {vec[0]: .6f}  {vec[1]: .6f}  {vec[2]: .6f}")