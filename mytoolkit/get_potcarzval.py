#!/usr/bin/env python
import os 
import re

from pprint import pprint
from collections import defaultdict

import numpy as np


# 读取POTCAR
with open('POTCAR', 'r') as f:
    lines = f.readlines()


# 获得元素名称
eleinfo = os.popen('grep -n "VRHFIN" POTCAR').readlines()
elements_names = []
for element in eleinfo:
    n = re.split(r'[:,\s,=]+', element)[2]
    elements_names.append(n)
# print(elements_names)

# 获得价电子个数
zvalinfo = os.popen('grep -B 1 "parameters from PSCTR" POTCAR').read()
zvalences = re.findall(r'\d+\.\d+', zvalinfo)
zvalences = [float(z) for z in zvalences]
# print(zvalences)

# 获得总电子轨道数以及文件中每个元素存储的起始行号
entries = os.popen('grep -n "entries" POTCAR').readlines()
eletronorbit_entries = []
eletron_orb_row   = []
for entry in entries:
    r = re.split(r'[:,\s]+', entry)[0]
    e_orb= re.split(r'[:,\s]+', entry)[1]
    eletron_orb_row.append(int(r))
    eletronorbit_entries.append(int(e_orb))
# print(eletronorbit_entries)
# print(eletron_orb_row)


# 获得每个元素的总电子轨道分布
eletron_configs = []
for r, e_orb in zip(eletron_orb_row, eletronorbit_entries):
    dstlines = lines[r+1:r+1+e_orb]
    new_dstlines = []
    for dl in dstlines:
        dl = dl.strip('\n').split()
        dl = [float(l) for l in dl]
        new_dstlines.append(dl)
    eletron_configs.append(new_dstlines[::-1])
# pprint(eletron_configs)

# 获得每个元素的价电子轨道分布
zval_eleconfigs = []
for eleconfig, elename, zval in zip(eletron_configs, elements_names, zvalences):
    calculated_zval = 0
    zval_eleconfig = []
    for idx, elecfg in enumerate(eleconfig):
        if calculated_zval < zval and not np.isclose(elecfg[-1], 0, rtol=1e-3):
            calculated_zval += elecfg[-1]
            zval_eleconfig.append(elecfg)
    zval_eleconfigs.append(zval_eleconfig)
# pprint(zval_eleconfigs)

# 获取轨道信息和轨道名称
angular_quantum = ['s', 'p', 'd', 'f']
element_eleorbit = defaultdict(dict)
for eleconfig, elename, zval in zip(zval_eleconfigs, elements_names, zvalences):
    zvalorbit_name_elenum = ''
    zvalorbit_name = ''
    eleconfig = sorted(eleconfig, key=lambda item: (item[0], item[1]))
    for elecfg in eleconfig:
        name_elenum = str(int(elecfg[0])) + angular_quantum[int(elecfg[1])] + str(elecfg[-1]) + ' '
        name = str(int(elecfg[0])) + angular_quantum[int(elecfg[1])] + " "
        zvalorbit_name_elenum += name_elenum
        zvalorbit_name += name

    element_eleorbit[elename]['orbitname'] = zvalorbit_name
    element_eleorbit[elename]['orbitelenum'] = zvalorbit_name_elenum
    element_eleorbit[elename]['zvalnumber'] = zval
    
pprint(dict(element_eleorbit))

# 打印信息
for ele_name, infos in element_eleorbit.items():
    print(ele_name)
    for infoname, info in infos.items():
        print("    {}".format(info))
print("For lobster input file basisfunctions")
for ele_name, infos in element_eleorbit.items():
    for infoname, info in infos.items():
        if infoname == 'orbitname':
            print("basisfunctions {:<4} {:<4}".format(ele_name, info))