#!/usr/bin/env python
import sys

from ase.io import read
import numpy as np

file = sys[1]
# 从文件中读取晶体结构
atoms = read(file)

# 计算距离矩阵
distances = atoms.get_all_distances()

# 将对角线元素设为1000
np.fill_diagonal(distances, 1000)

# 打印距离矩阵
print("最小距离是", np.min(distances))