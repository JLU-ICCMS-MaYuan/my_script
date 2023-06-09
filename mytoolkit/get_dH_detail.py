#!/usr/bin/env python3

import os
import re
import sys

print("Note: --------------------")
print("    这个脚本是从POTCAR中读取焓值 dH, 并且读取原子数 N, 然后, 计算每原子的焓值 dH/atom")
print("    所以你一定要确保你的OUTCAR是没有问题的")
print("    如果你想批量使用该脚本可以用这个命令: \n    for i in *; do if [ -d $i ]; then cd $i; echo $i; dHperatom.py ;cd ..; fi; done\n")

try:
    outcar_file = sys.argv[1]
except:
    outcar_file = "OUTCAR"

try:
    dH = os.popen(f"grep enthalpy {outcar_file}" + " | tail -n 1 | awk '{print $ 5}'").read().strip('\n')
    
    print("    dH = {:<12.8f} eV".format(float(dH)))
except:
    print("   OUTCAR有问题读不出来焓值")
    dH = 100000000000000.0
    N  = 1

begin_id = os.popen(f'grep -n "position of ions in cartesian coordinates" {outcar_file}').read().split(":")[0]
N = 0; row_id=int(begin_id)
while True:
    row_id = row_id+1
    content  = os.popen("sed -n '{}p' {}".format(row_id, outcar_file)).read().strip('\n')
    corrds   = re.findall(r"[-+]?\d+\.\d+", content)
    if len(corrds) == 3:
        N += 1
    else:
        break
print("    NumberofTotalAtoms = {:<12}".format(N))

print("    dH = {:<12.8f} eV/atom = {:<12.8f} meV/atom\n".format(float(dH)/N, float(dH)*1000/N))
