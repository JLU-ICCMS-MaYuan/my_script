#!/usr/bin/env python3

import sys
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import pandas as pd
import numpy as np
import collections

from scipy import interpolate

Enthalpy_curve_file = sys.argv[1]

Enthalpy_curve_data = pd.read_csv(Enthalpy_curve_file, index_col=0, header=0,)

new_Enthalpy_curve_data = Enthalpy_curve_data.dropna()

info="""创建一个目录叫:
D.dat
在该目录中写下要删除的内容
执行delete_pandas_columns.py formed_enthalpy.csv
生成new_formed_enthalpy.csv
"""
print(info)

with open("D.dat", 'r') as delete_dat:
    lines = delete_dat.readlines()


for line in lines:
# 获得要删除的列包含的关键词
    line = line.strip('\n')
    if line:
        # Enthalpy_curve_data = Enthalpy_curve_data.loc[:, ~Enthalpy_curve_data.columns.str.contains(line)]
        new_columns = []
        for col in new_Enthalpy_curve_data.columns:
            if line not in col: # 如果不包含 目标关键字, 就保留下来到new_columns
                new_columns.append(col)
            else:
                # print("delete {}".format(col))
                pass
        new_Enthalpy_curve_data = new_Enthalpy_curve_data[new_columns]
    else:
        break

new_Enthalpy_curve_data.to_csv("new_formed_enthalpy.csv")
print("完成")
# 2\*(6-Y1H6)+ 1\*(14-Sr1H22)