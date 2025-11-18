#!/usr/bin/env python3

from dpdata import System, LabeledSystem, MultiSystems

total_datas = MultiSystems()
# 53 76 77 
numlist = [i for i in range(1, 2)]
for i in numlist:
    print(f"----------convert outcar to raw and npy for VASP-{i}-----------")
    ls=LabeledSystem("../0.calypso-data/"+str(i)+"/"+"OUTCAR", fmt='outcar')
    total_datas.append(ls)

print(total_datas)
total_datas.to_deepmd_raw('1.trainset')
total_datas.to_deepmd_npy('1.trainset', set_size=5)
