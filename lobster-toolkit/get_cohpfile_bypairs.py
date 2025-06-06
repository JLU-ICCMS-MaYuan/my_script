#!/usr/bin/env python
import os
import sys
from pprint import pprint
from collections import defaultdict
from itertools import combinations, product

import numpy as np
from pymatgen.core.structure import Structure

info="""使用说明:
python get_pairs.py [COHPstartEnergy] [COHPendEnergy] [element1] [element2] [lower_d_limit] [upper_d_limit]
指定其实计算的能量COHPstartEnergy和终止计算的能量COHPendEnergy
指定两种元素element1和element2, 两种元素可以相同, 效果是求两种元素在单胞内所有原子间距离
指定原子间距离，获得指定距离下原子对。

获得指定COHPCAR.lobster中的费米面位置(第一列), COHP值(第二列), ICOHP值(第三列).
awk '$1 == "0.00000" {print $1, $2, $3}' COHPCAR.lobster

"""
print(info)

COHPstartEnergy = sys.argv[1]
COHPendEnergy   = sys.argv[2]
species_custom1 = sys.argv[3]
species_custom2 = sys.argv[4]

try:
    lower_d = float(sys.argv[5])
    upper_d = float(sys.argv[6])
except:
    lower_d = 0.0
    upper_d = 5.0


struct = Structure.from_file("POSCAR")

# 获得指定的两组元素的坐标
sites_custom1 = [site for idx, site in enumerate(struct.sites) if site.species_string == species_custom1] 
index_custom1 = [idx+1  for idx, site in enumerate(struct.sites) if site.species_string == species_custom1] # 因为python的索引编号是从0开始的
sites_custom2 = [site for idx, site in enumerate(struct.sites) if site.species_string == species_custom2] 
index_custom2 = [idx+1  for idx, site in enumerate(struct.sites) if site.species_string == species_custom2]

if species_custom1 == species_custom2:
    pairs = list(combinations(sites_custom1, r=2))
    idxs = list(combinations(index_custom1, r=2))
else:
    pairs = list(product(sites_custom1, sites_custom2))
    idxs = list(product(index_custom1, index_custom2))

# 获得每个原子对的距离并且存储
distances = []
for pair in pairs:
    d = struct.lattice.get_all_distances(pair[0].frac_coords, pair[1].frac_coords)[0][0]
    distances.append(np.round(d, decimals=3))

# 获得最终的数据
pairs_idxs_d = list(zip(pairs, idxs, distances))
pairs_idxs_d = sorted(pairs_idxs_d, key=lambda x: x[2])

# 打印相关信息
print("Note: --------------------")
print("{} {}".format("Number", "Elements"))
for idx, site in enumerate(struct.sites):
    print("{:<8} {:<8}".format(idx+1, site.specie))

print("Note: --------------------")
print("{} {} {} {}".format("Num", "Species", "Species", "distance"))
N_x = 0
for pair, idx, d in pairs_idxs_d:
    if lower_d <= d <= upper_d:
        N_x += 1
        print("{:<3}  {:>2}{:<4} {:>2}{:<4} {:>10.6f}".format(N_x , pair[0].specie, idx[0], pair[1].specie, idx[1], d))

with open("lobsterinIsmear_5", "w") as f:
    f.write('COHPstartEnergy  {}\n'.format(COHPstartEnergy)) 
    f.write('COHPendEnergy    {}\n'.format(COHPendEnergy))
    f.write('usebasisset pbeVaspFit2015\n') # 基组（直接使用根据vasp拟合的基组以及默认的基函数）
    for spe in struct.types_of_specie:
        f.write('basisfunctions {}\n'.format(spe.name)) # 基组（直接使用根据vasp拟合的基组以及默认的基函数）
    # 这里会出现一个非常严重的计算问题！！！！！！！！！！！！！！
    # lobster认为你指定的原子对是不具有周期性的
    # 你用pymatgen脚本找到的距离是包含周期性的，把这原子对输入给lobsterin
    # 它认不出来这个距离是周期性的，它会按照原胞内的距离考虑两个原子的成键。
    # 所以这里我抛弃了设置原子对来计算成键强度的方法。
    # 改用设置键长来获得原子对，lobster有自己的算法来获得原子对。
    for pair, idx, d in pairs_idxs_d:
        if lower_d <= d <= upper_d:
            f.write("cohpbetween atom {} and atom {}\n".format(idx[0], idx[1]))
    f.write("#cohpGenerator from {} to {} type {} type {}\n".format(lower_d, upper_d, species_custom1, species_custom2))

with open("lobsterinIsmear_0", "w") as f:
    f.write('COHPstartEnergy  {}\n'.format(COHPstartEnergy)) 
    f.write('COHPendEnergy    {}\n'.format(COHPendEnergy))
    f.write('usebasisset pbeVaspFit2015\n') # 基组（直接使用根据vasp拟合的基组以及默认的基函数）
    f.write('gaussianSmearingWidth 0.05\n')
    for spe in struct.types_of_specie:
        f.write('basisfunctions {}\n'.format(spe.name)) # 基组（直接使用根据vasp拟合的基组以及默认的基函数）
    # 这里会出现一个非常严重的计算问题！！！！！！！！！！！！！！
    # lobster认为你指定的原子对是不具有周期性的
    # 你用pymatgen脚本找到的距离是包含周期性的，把这原子对输入给lobsterin
    # 它认不出来这个距离是周期性的，它会按照原胞内的距离考虑两个原子的成键。
    # 所以这里我抛弃了设置原子对来计算成键强度的方法。
    # 改用设置键长来获得原子对，lobster有自己的算法来获得原子对。
    for pair, idx, d in pairs_idxs_d:
        if lower_d <= d <= upper_d:
            f.write("cohpbetween atom {} and atom {}\n".format(idx[0], idx[1]))
    f.write("#cohpGenerator from {} to {} type {} type {}\n".format(lower_d, upper_d, species_custom1, species_custom2))



dirs = "{}_{}_{}_{}".format(species_custom1, species_custom2, lower_d, upper_d)
if not os.path.exists(dirs):
    os.mkdir(dirs)
print("Note: ------------------------------------------------------")
print("    You need to do it before running Lobster-4.1.0")
print("    cp lobsterraw lobsterin")
print("    You need to do it after running Lobster-4.1.0")
print("    cp {COBICAR.lobster,COHPCAR.lobster,COOPCAR.lobster,ICOBILIST.lobster,ICOHPLIST.lobster,ICOOPLIST.lobster,lobsterin}" + "  {}".format(dirs))
print("------------------------------------------------------------")
