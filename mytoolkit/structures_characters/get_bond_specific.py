#!/usr/bin/env python
import os
import sys
from pprint import pprint
from collections import defaultdict


import numpy as np
from pymatgen.core.structure import Structure
from pymatgen.core.bonds import CovalentBond


info="""使用说明:
python get_bond.py [element1] [element2]
指定两种元素element1和element2, 两种元素可以相同
效果是求两种元素在单胞内所有原子间距离
"""
print(info)

species_custom1 = sys.argv[1]
species_custom2 = sys.argv[2]

struct = Structure.from_file("POSCAR")

# 获得指定的两组元素的坐标
frac_custom1 = [site.frac_coords for site in struct.sites if site.species_string == species_custom1] 
frac_custom2 = [site.frac_coords for site in struct.sites if site.species_string == species_custom2] 

# 获取指定的两组元素的距离, 奖距离矩阵转化为一个一维数组
bondslength = struct.lattice.get_all_distances(frac_custom1, frac_custom2)
if species_custom1 == species_custom2:
    # 如果求相同的两个元素的距离, 需要获得其对称矩阵并求上三角矩阵
    bondslength_rmdup = np.triu(bondslength, k=1).flatten() # k=0表示主对角线的位置，k=1表示主对角右移1，k=-1表示对角线左移1
    bondslength_rmdup = np.round(bondslength_rmdup[bondslength_rmdup > 0.001], decimals=3)
    countting_total_bonding_number = len(bondslength_rmdup)
    estimating_total_bonding_number = len(frac_custom1)*(len(frac_custom1)-1)/2
else:
    # 如果求不同相同的两个元素的距离, 只需去掉小于等于0的数即可
    bondslength_rmdup = np.round(bondslength[bondslength > 0.001], decimals=3)
    countting_total_bonding_number = len(bondslength_rmdup)
    estimating_total_bonding_number = len(frac_custom1)*len(frac_custom2)

# 统计每种长度的键的数目并存储在字典bondinglength_numbers中
bondinglength_numbers = defaultdict(int)
for bd1 in set(bondslength_rmdup):
    count = np.where(np.isclose(bd1, bondslength_rmdup, rtol=1e-4), 1, 0)
    # 这里判断两个键长是否相同的条件需要非常苛刻, 因为键长保留到小数点后3位, 
    # 所以两个键长的差别最小要通过小数点后第4位体现, 所以起码应该保留小数点后4位
    bondinglength_numbers[np.round(bd1, decimals=3)] = sum(count)
bondinglength_numbers = sorted(bondinglength_numbers.items())

# 打印相关信息
print("Note: --------------------")
print("Total number of bonding by estimate: {:<10}".format(int(estimating_total_bonding_number)))
print("Total number of bonding by countting: {:<10}".format(int(countting_total_bonding_number)))
print("{} different bonding length".format(len(set(bondslength_rmdup))))
print("Count the number of each type of bonding length")
print("{:>10} {:>10}".format("Length(A)", "Number"))
summing_total_bonding_number = 0
for bd, num in bondinglength_numbers:
    print("{:>10.3f} {:>10}".format(bd, num))
    summing_total_bonding_number += num
print("Total number of bonding by summing: {:<10}".format(int(summing_total_bonding_number)))
    