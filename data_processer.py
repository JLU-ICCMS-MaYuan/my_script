#!/usr/bin/env python

import pandas as pd
import os

# 强制pandas不使用第一列作为索引（行名） index_col=False
df1 = pd.read_csv('convexhull.csv', index_col=False) 

# 删除重复项
# keep：有三个可选参数，分别是 first、last、False，
# 默认为 first，表示只保留第一次出现的重复项，删除其余重复项，
# last 表示只保留最后一次出现的重复项，
# False 则表示删除所有重复项。
# inplace, 默认为 False 表示删除重复项后返回一个副本，若为 Ture 则表示直接在原数据上删除重复项。
df1.drop_duplicates(subset='Number', keep='first', inplace=True)

# 处理有NaN这样缺失值的行, inplace表示直接在原DataFrame修改 
# 如果有NaN, 那么axis=0表示删除存在NaN的一整行
# inplace, 默认为 False 表示删除重复项后返回一个副本，若为 Ture 则表示直接在原数据上删除重复项。
df1.dropna(axis=0, inplace=True) 

# 删除焓值为 610612509 的行
# index = df1[(df1.enthalpy - 610612509.0 < 0.1)].index.tolist() 表示返回与610612509差值小于0.1的所有行的行索引，将这些行索引存储为列表
df1.drop(
    index=df1[(abs(df1.enthalpy - 610612509.0) < 0.1)].index.tolist(),
    inplace=True
    )

df1.to_csv('nnconvexhull.csv', index=False)

