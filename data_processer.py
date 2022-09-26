#!/usr/bin/env python

import pandas as pd
import os

_tot_list = os.popen("grep 'NULL' dir_*/Analysis_Output.dat").read()
tot_list = _tot_list.replace('(', ' ')
tot_list = tot_list.replace(')', ' ')

NullNum = []
for temp in tot_list.split('\n'):
    info = temp.split()
    if len(info) != 0:
        NullNum.append(info[2])

NullNum = list(map(int, NullNum))

df1 = pd.read_csv('convexhull.csv', index_col=False) # 强制pandas不使用第一列作为索引（行名）
df1.drop_duplicates(subset='Number', keep='first', inplace=True)
df1.dropna(axis=0, inplace=True) # 处理有NaN这样缺失值的行, inplace表示直接在原DataFrame修改 
df1.to_csv('nnconvexhull.csv', index=False)

