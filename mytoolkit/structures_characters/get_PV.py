#!/usr/bin/env python3

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-p',  '--pressure', type=float, default=None, help='pressure with unit is GPa')
parser.add_argument('-v',  '--volume',   type=float, default=None, help='volume   with unit is angstrom^3')
parser.add_argument('-pv', '--pv-item',  type=float, default=None, help='PV iterm with unit is eV')

args = parser.parse_args()

v = args.volume
p = args.pressure
pv = args.pv_item
GPaA3_2_eV = 6.242e-3

# 计算PV项能量
if (pv is None) and (v is not None) and (p is not None):  
    PV_item = p * v * GPaA3_2_eV  # 计算 PV 项能量
    print(f"Calculated PV item (eV): {PV_item}")
    
# 计算压强
elif (p is None) and (v is not None) and (pv is not None):
    p = pv / (v * GPaA3_2_eV)  # 根据 PV 项和体积计算压强
    print(f"Calculated pressure (GPa): {p}")

# 计算体积
elif (v is None) and (p is not None) and (pv is not None):
    v = pv / (p * GPaA3_2_eV)  # 根据 PV 项和压强计算体积
    print(f"Calculated volume (Å³): {v}")

# 如果都有或者都没有输入，打印输入参数
else:
    print(f"Invalid or incomplete input. Please provide two values:")
    print(f"Pressure (p): {p}, Volume (v): {v}, PV item (pv): {pv}")