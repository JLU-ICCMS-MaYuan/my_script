#!/usr/bin/env python3

from dpdata import System, LabeledSystem, MultiSystems

total_datas = MultiSystems()
# 53 76 77 
numlist = [i for i in range(1, 501)]
for i in numlist:
    print(f"----------convert outcar to raw and npy for VASP-{i}-----------")
    ls=LabeledSystem(str(i)+"/"+"OUTCAR", fmt='outcar')
    total_datas.append(ls)

total_datas.to_deepmd_raw('../dp-rawdata')
total_datas.to_deepmd_npy('../dp-rawdata', set_size=5)
