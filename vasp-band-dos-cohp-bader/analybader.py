#!/usr/bin/env python3

import os
import sys
import pandas as pd
from numpy import round
from pprint import pprint

print("You have to gaurantee there are bader chgsum.pl and analybader")

print("Note: --------------------")
print("Step 1, you will get CHGCAR_sum")
print("./chgsum.pl AECCAR0 AECCAR2")
os.system("./chgsum.pl AECCAR0 AECCAR2")

print("Note: --------------------")
print("Step 2, you will get ACF.dat AVF.dat BCF.dat, ACF.dat ")
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
info="""please input parameters by the order [elementsname], [numberofatoms] [zvalforelement], such as,
    python analybader.py La 4 11 H 40 1
the order of elements must be consistent with that of POTCAR"""

ele_natom_zval = sys.argv[1:]
if not ele_natom_zval:
    print(info)
    print("please input right parameters")
    sys.exit(1)
else:
    print(info)
    print(ele_natom_zval)

ele_natom_zval = [
        {
        "name":ele_natom_zval[y], "natom":int(ele_natom_zval[y+1]), "zval":int(ele_natom_zval[y+2])
        } for y in range(0, len(ele_natom_zval), 3)
        ]
pprint(ele_natom_zval)

serialnumber = 0
ele_serialnum_zval = []
for EleAtmZ in ele_natom_zval:
    for i in range(1, EleAtmZ['natom']+1, 1):
        ele_serialnum_zval.append([i+serialnumber, EleAtmZ['name'], EleAtmZ['name']+str(i), EleAtmZ['zval']])
    else:
        serialnumber += i

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
