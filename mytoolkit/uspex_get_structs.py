#!/usr/bin/env python3

import os
import sys

print("Note:--------------------")
print("    该脚本的作用是从uspex中提取出指定fitness能量范围内的结构")
print("    一旦指定了fitness的范围, 脚本就会从extended_convex_hull中获得相应能量范围内的结构的编号")
print("    然后根据指定的编号在symmetrized_structures.cif文件中获得相应的cif结构")
print("    最后将获得的结构存储在seeds.vasp中")

lower_limit = float(sys.argv[1])
upper_limit = float(sys.argv[2])

with open("extended_convex_hull", "r") as ech:
    content1 = ech.readlines()
content1 = content1[6:]

struct_ids = []
for line in content1:
    line = line.strip("\n").split()
    fitness = float(line[8])
    if fitness<=upper_limit and fitness>=lower_limit:
        struct_id = line[0]
        struct_ids.append(struct_id)

with open("symmetrized_structures.cif", "r") as cif:
    content2 = cif.readlines()

for line in content2:
    input(line.strip('\n'))
