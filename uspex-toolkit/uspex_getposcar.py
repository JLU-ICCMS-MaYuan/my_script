#!/usr/bin/env python3

import os
import sys
info = '''Note:--------------------
    该脚本的作用是从uspex中提取出指定fitness能量范围内的结构
    一旦指定了fitness的范围, 脚本就会从extended_convex_hull中获得相应能量范围内的结构的编号
    指定的fitness范围单位是eV/atom, 所以如果指定范围是0~0.010eV, 那么这么写: uspex_getposcar.py 0 0.01
    获得的结构将存储在poscars目录中, 更加详细的写法有：
uspex_getposcar.py 0 0.1 LaScH_ True  True  True
uspex_getposcar.py 0 0.1 LaSc_  True  True  False
uspex_getposcar.py 0 0.1 LaH_   True  False True 
uspex_getposcar.py 0 0.1 ScH_   False True  True
uspex_getposcar.py 0 0.1 La_    True  False False
uspex_getposcar.py 0 0.1 Sc_    False True  False
uspex_getposcar.py 0 0.1 H_     False False True
'''

print(info)


lower_limit = float(sys.argv[1])
upper_limit = float(sys.argv[2])
try:
    prefix = sys.argv[3]
    first_ele_include = eval(sys.argv[4])
    second_ele_include= eval(sys.argv[5])
    third_ele_include = eval(sys.argv[6])
    print('try')
    print(first_ele_include, second_ele_include, third_ele_include)
    print(sys.argv[:])
except:
    prefix = ''
    first_ele_include  = True
    second_ele_include = True 
    third_ele_include  = True 
    print(first_ele_include, second_ele_include, third_ele_include)


# 从extended_convex_hull中获得相应结构的编号
with open("extended_convex_hull", "r") as ech:
    content1 = ech.readlines()
content1 = content1[6:]

struct_ids = []
for line in content1:
    line = line.strip("\n").split()
    fitness = float(line[8])
    composition = list(map(int, line[2:5]))
    if fitness<=upper_limit and fitness>=lower_limit:
        if  (first_ele_include  and composition[0]>0) and \
            (second_ele_include and composition[1]>0) and \
            (third_ele_include  and composition[2]>0):
            struct_id = line[0]
            struct_ids.append(struct_id)
        elif (not first_ele_include and composition[0]==0) and \
             (second_ele_include and composition[1]>0) and \
             (third_ele_include  and composition[2]>0):
            struct_id = line[0]
            struct_ids.append(struct_id)
        elif (first_ele_include  and composition[0]>0) and \
             (not second_ele_include and composition[1]==0) and \
             (third_ele_include  and composition[2]>0):
            struct_id = line[0]
            struct_ids.append(struct_id)
        elif (first_ele_include  and composition[0]>0) and \
             (second_ele_include and composition[1]>0) and \
             (not third_ele_include and composition[2]==0):
            struct_id = line[0]
            struct_ids.append(struct_id)
        elif (first_ele_include and composition[0]>0) and \
             (not second_ele_include and composition[1]==0) and \
             (not third_ele_include  and composition[2]==0):
            struct_id = line[0]
            struct_ids.append(struct_id)
        elif (not first_ele_include  and composition[0]==0) and \
             (second_ele_include and composition[1]>0) and \
             (not third_ele_include  and composition[2]==0):
            struct_id = line[0]
            struct_ids.append(struct_id)
        elif (not first_ele_include  and composition[0]==0) and \
             (not second_ele_include and composition[1]==0) and \
             (third_ele_include  and composition[2]>0):
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