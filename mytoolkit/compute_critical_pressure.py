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
# print(Enthalpy_curve_data.head(5))
# 去除NaN值
Enthalpy_curve_data = Enthalpy_curve_data.dropna()

critical_press_data = pd.DataFrame()
for chemical_path, press_enthalpy in Enthalpy_curve_data.iteritems():
    # print(f"chemical_path={chemical_path}")
    press_data = np.array(press_enthalpy.index)
    enth_data  = np.array(press_enthalpy.values)

    # 定义一个函数，返回曲线数据的差值
    tck = interpolate.make_interp_spline(x=press_data, y=enth_data, k=1)

    # 使用root_scalar函数来寻找根
    piecewise_polynomial = interpolate.PPoly.from_spline(tck, extrapolate=None)
    critical_press = piecewise_polynomial.roots()
    critical_press = critical_press[np.where(np.logical_and(critical_press>=-100, critical_press<=300))]

    if len(critical_press) == 1:
        critical_press_data.at[chemical_path, "lower_limit"] = critical_press[0]
        critical_press_data.at[chemical_path, "upper_limit"] = 100000
    elif len(critical_press) == 2:
        critical_press_data.at[chemical_path, "lower_limit"] = critical_press[0]
        critical_press_data.at[chemical_path, "upper_limit"] = critical_press[1]
    
    # critical_press_data = critical_press_data.sort_values(by=["lower_limit", "upper_limit"])
# 获取零点的值

print(critical_press_data)

# print(f"chemical_path={critical_press_data['lower_limit'].nlargest(4).index[:]} \
#         lower_limit={critical_press_data['lower_limit'].nlargest(4).iloc[:]}\n")
print(f"lower_limit=\n{critical_press_data['lower_limit'].nlargest(12).iloc[:]}\n")
# print(f"chemical_path={critical_press_data['upper_limit'].nsmallest(4).index[:]} \
#         upper_limit={critical_press_data['upper_limit'].nsmallest(4).iloc[:]}\n")
# 13*(1-H1)+ 1*(3-La1H3)+ 1*(4-Y1H2)+ 1*(9-Ce1H3)+ 1*(13-Th1H3)+ 4*(19-Be1H2) -> 1*(20-La1Ce1Y1Th1Be4H32)