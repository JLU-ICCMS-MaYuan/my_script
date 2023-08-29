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

Enthalpy_curve_data = Enthalpy_curve_data.loc[:, ~Enthalpy_curve_data.columns.str.contains(name)]


Enthalpy_curve_data.to_csv("new-formed-enthalpy.csv")
