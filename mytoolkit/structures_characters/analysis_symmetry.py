#!/usr/bin/env python3
import sys
import time
from pathlib import Path

import pandas as pd
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer



try:
    print("使用格式:  symbolspg.py 目录路径 对称性的精度 指定你要分析的空间群号")
    print("例如:  symbolspg.py ./ 0.01 229")
    path = Path(sys.argv[1])
    symtol = eval(sys.argv[2])
    specified_spgnum = eval(sys.argv[3])
except:
    print("如果指定路径和寻找对称性精度时出现问题就会需要手动重新输入")
    path   = Path(input("请输入待计算的目录路径:\n"))
    symtol = eval(input("请输入寻找对称性精度:\n"))
    specified_spgnum = eval(input("指定你要分析的空间群号:\n"))
print("对称性的精度: {}".format(str(symtol)))
print("批量寻找的结构的名称为POSCAR_*，请注意是否需要修改")
struct_name = input("请输入结构名称的正则表达式, 默认是POSCAR_*:\n")
if not struct_name:
    struct_name = "POSCAR_*"
print("批量寻找的结构的名称为: {}".format(struct_name))
time.sleep(3)

with open(str(symtol)+"_analysis_result", "w") as f:
    f.write("{:<30}   {:<8}   {:<8}   {:<10}\n".format("file_name", "spg_num", "spg_symbol", "formula"))
    for file in path.glob(struct_name):
        file_name = file.name
        print("当前处理的文件为: {}".format(file_name))
        struct = Structure.from_file(file)
        name = struct.composition.get_integer_formula_and_factor()[0]
        spg_num = SpacegroupAnalyzer(struct, symprec=symtol).get_space_group_number()
        spg_symbol = SpacegroupAnalyzer(struct).get_space_group_symbol()
        # print("{:<30}   {:<5}   {:<8}   {:<10}".format(file_name, spg_num, spg_symbol, name))
        f.write("{:<30}   {:<5}   {:<8}   {:<10}\n".format(file_name, spg_num, spg_symbol, name))


data = pd.read_table(str(symtol)+"_analysis_result", sep="\s+")

# 计算 spg_num 列中每个值的数量
spg_counts = data['spg_num'].value_counts()

# 计算 总的结构数
total_number  = len(data)

# 计算 指定要分析的空间群号 的数量
try:
    specified_count = spg_counts.loc[specified_spgnum]
    poscar_names = data.loc[data['spg_num'] == specified_spgnum, 'file_name']
    print("指定的空间群的数目{}/{}, {:<4.3f}".format(specified_count, total_number, specified_count/total_number*100))
    print("Their names are")
    print(poscar_names)
except Exception as e:
    print("报错, 也许是这些结构中没有这个空间群对称性")
    print(e)

