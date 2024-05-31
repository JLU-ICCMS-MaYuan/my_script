#!/usr/bin/env python3
import sys

import pandas as pd

'''
csv_merge.py 二元数据.csv 三元数据.csv 四元数据.csv 五元数据.csv 六元数据.csv
'''

filenames = sys.argv[1:]

all_dfs = []
for fn in filenames:
    print(fn)
    one_df = pd.read_csv(fn, index_col=0, header=0, na_values=['--'])
    all_dfs.append(one_df)

# 使用 concat 函数将它们 axis=0 按列合并 
# axis = 0指的是拼接行（向下拼接）
# axis = 1指的是拼接列（向右拼接）

# 如果 ignore_index 设置为 True，
# 则合并后的数据帧将重新生成一个新的索引（默认为从 0 开始的整数索引），
# 忽略原始数据帧的索引。这对于合并多个数据帧，
# 并将它们整合成一个新的数据帧，不考虑原始索引的情况很有用。

# 如果 ignore_index 设置为 False（默认值），
# 则合并后的数据帧将保留原始数据帧的索引。
# 这意味着合并后的数据帧将保留每个原始数据帧的索引，
# 并在合并后的数据帧中可能存在重复的索引值。
merged_df = pd.concat(all_dfs, axis=1, ignore_index=False)
merged_df.to_csv('origin.csv', header=True, index=True)
