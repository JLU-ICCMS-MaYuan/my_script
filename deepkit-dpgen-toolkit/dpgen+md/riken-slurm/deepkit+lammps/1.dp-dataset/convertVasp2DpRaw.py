#!/usr/bin/env python3

from dpdata import System, LabeledSystem, MultiSystems

total_datas = MultiSystems()
# 53 76 77 
numlist = [
'../2000K/1.potim_0.1_and_smass_2',
'../2000K/2.potim_0.1_and_smass_0',
'../2000K/3.potim_0.1_and_smass_0',
'../2500',
'../2500/2.potim_0.1_and_smass_0',
'../3000K/1.potim_0.1_and_smass_2',
'../3000K/2.potim_0.1_and_smass_0',
'../3500',
'../4000K/1.potim_0.1_and_smass_2',
'../4000K/2.potim_0.1_and_smass_0',
]
for i in numlist:
    print(f"----------convert outcar to raw and npy for VASP-{i}-----------")
    try:
        ls=LabeledSystem(str(i)+"/"+"OUTCAR", fmt='outcar')
        total_datas.append(ls)
    except ValueError:
        print(f'error for {i}')

total_datas.to_deepmd_raw('1.traindata')
total_datas.to_deepmd_npy('1.traindata', set_size=5)
