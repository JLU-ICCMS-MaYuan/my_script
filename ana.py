#!/usr/bin/env python
import sys
import pandas as pd
import matplotlib.pyplot as plt


path = sys.argv[1]
# 读取文件并跳过前3行
data = pd.read_table(path, skiprows=3, header=None, names=["struct_name", "spg_num", "spg_symbol", "formula"], sep="\s+")

# 计算 spg_num 列中每个值的数量
spg_counts = data['spg_num'].value_counts()


# 计算 spg_num 为 200 的比例
total_number  = len(data)
dst_number = spg_counts.loc[221]
dst_ratio = dst_number / total_number
rest_number   = total_number - dst_number
print(dst_number)

# plt.bar(1, dst_number, label='221',fc = 'y')
# plt.bar(1, rest_number, label='other',fc = 'r')
# plt.show()