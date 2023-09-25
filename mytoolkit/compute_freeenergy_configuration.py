#!/usr/bin/env python3

import sys

import pandas as pd

# 第一步：读入数据formed-enthalpy.csv
# 第二步：给每一列都加上TS项
# 第三步；输出formed-freeenergy.csv

if __name__ == "__main__":
    
    Enthalpy_curve_file = sys.argv[1] # 单位是meV/atom
    T_times_entropy     = float(sys.argv[2]) # 单位是meV/atom
    Enthalpy_curve_data = pd.read_csv(Enthalpy_curve_file, index_col=0, header=0,)

    free_energy_data    = Enthalpy_curve_data + T_times_entropy

    # print(result_dt)
    total_result_pd = pd.DataFrame(data=free_energy_data,) #index=0, columns=0)
    total_result_pd.to_csv("formed_freeenergy.csv")
    



