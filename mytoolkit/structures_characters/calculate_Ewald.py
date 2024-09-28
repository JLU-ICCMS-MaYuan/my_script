#!/usr/bin/env python3

import os
import re
import argparse  

from pymatgen.core import Structure
from pymatgen.analysis.energy_models import EwaldSummation
from pymatgen.io.vasp import Poscar

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

# 添加 -v 参数，使用自定义的action来解析  
parser.add_argument('-v', '--valences', type=str, nargs='+', required=True, help='Valence pairs in format "Element +Valence"') 

# 解析命令行参数  
args = parser.parse_args()  

valences = parse_valences(args.valences)

# 获取当前路径下所有POSCAR_开头的文件
poscar_files = [f for f in os.listdir('.') if re.match(r'POSCAR_.*', f)]

# 存储文件名和对应的静电能量
energy_data = []

for poscar_file in poscar_files:
    # 读取结构
    structure = Poscar.from_file(poscar_file).structure
    # 计算静电能量
    structure.add_oxidation_state_by_element(valences)
    energy = EwaldSummation(structure).total_energy
    # 记录文件名和能量
    energy_data.append((poscar_file, energy))

# 按能量从低到高排序
energy_data.sort(key=lambda x: x[1])

# 将结果写入electrostatic_energy.dat文件
with open('electrostatic_energy.dat', 'w') as f:
    for poscar_name, energy in energy_data:
        f.write(f"{poscar_name} {energy}\n")

