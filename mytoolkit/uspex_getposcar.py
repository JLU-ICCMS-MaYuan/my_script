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
try:
    prefix = sys.argv[3]
except:
    prefix = ''


# 从extended_convex_hull中获得相应结构的编号
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
titlenames = ["EA"+idx for idx in struct_ids]
print("Note:--------------------")
print("    在extended_convex_hull中的fitness符合你的能量要求的结构有{}个".format(len(titlenames)))
# print("    它们的编号分别是：")
# for idx in struct_ids:
#     print(idx)

# 从相应结构的编号对应获得symmetrized_structures中获得相应结构的title的行号
with open("gatheredPOSCARS", "r") as poscars:
    content2 = poscars.readlines()
struc_ids = []
for lineid, line in enumerate(content2):
    for title in titlenames:
        if title+' ' in line:
            struc_ids.append([lineid, line.strip('\n')])
print("Note:--------------------")
print("    在gatheredPOSCARS找到的符合你的能量要求的结构有{}个".format(len(struc_ids)))
# print("    它们的编号分别是：")
# for title in struc_ids:
#     print(title)
# print("    如果在gatheredPOSCARS找到的符合你的能量要求的结构数 少于 在extended_convex_hull中的fitness符合你的能量要求的结构数，可能是因为你的续算导致部分结构在其它results*中")


# 根据相应结构的行号获得所有的结构，并将它们都存储在ciffiles目录中
poscar_s = []
for lineid, line in struc_ids:
    poscar = []
    for structinfo in content2[lineid:]:
        if "EA" in structinfo and structinfo.strip('\n') != line:
            break
        else:
            poscar.append(structinfo)
    poscar_s.append(poscar)
print("Note:--------------------")
print("    在{}提取出{}个poscar结构".format("gatheredPOSCARS", len(poscar_s)))
if os.path.exists("poscars"):
    os.system("rm -fr poscars")
    os.system("mkdir poscars")
else:
    os.system("mkdir poscars")
for poscar in poscar_s:
    with open(os.path.join("poscars", prefix+poscar[0].strip('\n').split()[0]+".vasp"), "w") as f:
        f.writelines(poscar)
