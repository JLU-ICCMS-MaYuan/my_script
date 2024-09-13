#!/usr/bin/env python3
import os
from dpdata import System, LabeledSystem, MultiSystems

print("You had better run it in `1.dp-data/2.testset`")
current_directory = os.getcwd()
subdirectories = [d for d in os.listdir(current_directory) if os.path.isdir(os.path.join(current_directory, d))]
total_datas = MultiSystems()
# 遍历所有子目录  
totalnum = 0
for subdir in subdirectories:
    ls=LabeledSystem(subdir,  fmt="deepmd/npy")
    totalnum+=len(ls)
    total_datas.append(ls)
print(totalnum)
