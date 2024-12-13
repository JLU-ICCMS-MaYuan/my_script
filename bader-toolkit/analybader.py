#!/usr/bin/env python3

# 学习脚本
# https://zhuanlan.zhihu.com/p/541611145
import os
import argparse

import pandas as pd
from numpy import round
from pprint import pprint

# 设定 argparse 参数
parser = argparse.ArgumentParser(description="Bader charge analysis script. For example:  python analybader.py -e La Th Ce Be H -n 1 1 1 3 24 -z 11 12 12 4 1")
parser.add_argument("-e", "--elements",       type=str, nargs='+', help="Element names")
parser.add_argument("-n", "--numberofatoms",  type=int, nargs='+', help="Number of atoms")
parser.add_argument("-z", "--zvalence",       type=int, nargs='+', help="Zvalence.")
args = parser.parse_args()

print("You have to gaurantee there are bader chgsum.pl and analybader")

print("Note: --------------------")
print("Step 1, you will get CHGCAR_sum")
print("./chgsum.pl AECCAR0 AECCAR2")
os.system("./chgsum.pl AECCAR0 AECCAR2")

print("Note: --------------------")
print("Step 2, you will get ACF.dat AVF.dat BCF.dat")
print("./bader CHGCAR -ref CHGCAR_sum")
os.system("./bader CHGCAR -ref CHGCAR_sum")

print("Note: --------------------")
print("Step 3, The program will tell you general situation")
charge_list = []
with open('ACF.dat', 'r') as infile:
    lines = infile.readlines()
    for line in lines:
        if line.rstrip().split()[0].isdigit():
            charge_list.append(line.split())


ele_natom_zval = [
        {
        "name":name, "natom":natom, "zval":zval
        } for name, natom, zval in zip(args.elements, args.numberofatoms, args.zvalence)
        ]
pprint(ele_natom_zval)

serialnumber = 0
ele_serialnum_zval = []
for EleAtmZ in ele_natom_zval:
    for i in range(1, EleAtmZ['natom']+1, 1):
        ele_serialnum_zval.append([i+serialnumber, EleAtmZ['name'], EleAtmZ['name']+str(i), EleAtmZ['zval']])
    else:
        serialnumber += i
pprint(ele_serialnum_zval)

new_ele_tranferE = []
for charge, EleAtmZ in zip(charge_list, ele_serialnum_zval):
    if int(charge[0]) == int(EleAtmZ[0]):
        transferE = round(float(charge[4])-float(EleAtmZ[3]), decimals=6)
        new_ele_tranferE.append([EleAtmZ[1], EleAtmZ[2], transferE])

df = pd.DataFrame(new_ele_tranferE, columns=['elements', 'labels', 'transferE'])
classify = df.groupby('transferE')['labels'].apply(list)
print("To classify labels in according with transferE")
print(classify)


meanvalue = df.groupby('elements')['transferE'].mean().round(6)
print("The average gain and loss of charge for each element")
pprint(meanvalue)
