#!/work/home/mayuan/miniconda3/envs/cage/bin/python3

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

df2 = pd.DataFrame({'Number': NullNum})
df2.drop_duplicates(subset='Number', keep='first', inplace=True)
print(df2)
df1 = pd.read_csv('convexhull.csv', index_col=False)
df1.drop_duplicates(subset='Number', keep='first', inplace=True)
print(df1)
df3 = df1.append(df2)
df3.drop_duplicates(subset='Number', keep=False, inplace=True)
df3.dropna(axis=0, inplace=True)
print(df3)
df3.to_csv('nnconvexhull.csv', index=False)

