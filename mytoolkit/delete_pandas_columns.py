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

Enthalpy_curve_data = Enthalpy_curve_data.dropna()

# 获得要删除的列包含的关键词
name = input("获得要删除的列包含的关键词\n")


# Enthalpy_curve_data = Enthalpy_curve_data.loc[:, ~Enthalpy_curve_data.columns.str.contains(name)]
new_colums = []
for col in Enthalpy_curve_data.columns:
    if name not in col:
        new_colums.append(col)
    else:
        # print("delete {}".format(col))
        pass

# print(new_colums);input()
new_Enthalpy_curve_data = Enthalpy_curve_data[new_colums]

new_Enthalpy_curve_data.to_csv("new-formed-enthalpy.csv")

# 2\*(6-Y1H6)+ 1\*(14-Sr1H22)