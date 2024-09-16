#!/usr/bin/env python3

import sys
import argparse  
import re  

from pymatgen.core import Structure, Lattice
from pymatgen.analysis.energy_models import EwaldSummation
from pymatgen.transformations.standard_transformations import SupercellTransformation
from pymatgen.io.vasp import Poscar

import itertools
import numpy as np

from multiprocessing import Pool, cpu_count

def parse_valences(arg) -> dict:  
    """  
    解析化合价参数，格式为 "元素 +化合价" 对，返回字典。  
    """  
    valences = {}  
    # print(arg)
    arg = [arg[i:i + 2] for i in range(0, len(arg), 2)]
    # print(arg)
    for item in arg:  
        valences[item[0]] = float(item[1])  
    return valences  

# 创建解析器  
parser = argparse.ArgumentParser(description="Process input parameters.")  
# 添加 -i 参数  
parser.add_argument('-i', '--inputfilename', type=str, required=True, help='Input filename')  
# 添加 -v 参数，使用自定义的action来解析  
parser.add_argument('-v', '--valences', type=str, nargs='+', required=True, help='Valence pairs in format "Element +Valence"') 
# 添加 -d 参数，现在它应该接收三个整数  
parser.add_argument('-d', '--dimforsupercell', nargs=3, type=int, default=[1,1,1], help='Three integers representing the dimensions for the supercell (x, y, z)')  
# 添加 -r 参数，用于指定删除的氢原子的数量  
parser.add_argument('-r', '--remove_hydrogens', type=int, default=0, help='Number of hydrogen atoms to remove (default: 0)')  
# 添加 -n 参数，用于指定输出静电能最低的n个结构  
parser.add_argument('-n', '--num_structures', type=int, default=10, help='Number of structures with the lowest electrostatic energy to output (default: 1)')

# 解析命令行参数  
args = parser.parse_args()  
inputfilename = args.inputfilename
valences = parse_valences(args.valences)
dimforsupercell = args.dimforsupercell
remove_hydrogens = args.remove_hydrogens
num_structures = args.num_structures
info = '''what you input:
inputfilename: {}
    valences: {}
    dimforsupercell: {}
    remove_hydrogens: {}
    num_structures: {}

You can use it by:
    delete_Hatoms.py -i ../../LaSc2H24.vasp -v La 3 Sc 3 H -0.375 -d 1 1 1 -r 1 -n 50
'''.format(inputfilename, valences, dimforsupercell, remove_hydrogens, num_structures)
print(info)

s1 = Structure.from_file(inputfilename)

supercell = sys.argv[2:5]
# sc = SupercellTransformation().from_scaling_factors(dimforsupercell[0], dimforsupercell[1], dimforsupercell[2])
sc = SupercellTransformation().from_scaling_factors(*dimforsupercell)
s2 = sc.apply_transformation(s1)
# print(s2)

# 找到所有氢原子的位置
h_indices = [i for i, site in enumerate(s2) if site.species_string == "H"]

# 计算需要移除的氢原子数量
print('delete {} hydrogens from {} hydrogens'.format(remove_hydrogens, len(h_indices)))

# 生成所有可能的排列组合，删除氢原子
combinations = list(itertools.combinations(h_indices, remove_hydrogens))
print("Combinational situations (r={}): {}".format(remove_hydrogens, len(combinations)))

# 定义一个列表来存储结果
results  = []
Lowest_structures = []

# 遍历每一种组合，计算静电能
for combination in combinations:
    s3 = s2.copy()
    s3.remove_sites(combination)
    s3.add_oxidation_state_by_element(valences)
    energy = EwaldSummation(s3).total_energy
    #print(energy)
    results.append((combination, energy))
    if len(Lowest_structures) <= num_structures:
        Lowest_structures.append([s3, energy])
        Poscar(s3).write_file('POSCAR_'+str(len(Lowest_structures)))
    else:
        max_energy = max(Lowest_structures, key=lambda x: x[1])[1]
        if energy < max_energy:
            for i in range(len(Lowest_structures)):
                if np.isclose(Lowest_structures[i][1], max_energy):
                    Lowest_structures[i] = [s3, energy]
                    Poscar(s3).write_file('POSCAR_'+str(i+1))
                    break

Lowest_structures.sort(key=lambda x: x[1])

with open('electrostatic_energy.dat', 'w') as f:
    for idx, (s, energy) in enumerate(Lowest_structures):
        print('{:<}  {:>10.6f}'.format(idx+1, energy), file=f)
        Poscar(s).write_file('POSCAR_'+str(idx+1))