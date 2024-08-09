#!/usr/bin/env python3
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# 从命令行参数读取文件名
filename = "lambda.out"

# 读取文件内容
with open(filename, 'r') as file:
    data = file.readlines()

# 常量定义：1 Ry = 13.605693 eV
Ry2ev = 13.605693

# 数据处理
degauss = []
N_Ef_Ry = []
N_Ef_eV = []

# 忽略 Tc 行及其以下的内容
skip = False
for line in data:
    if 'T_c' in line:
        skip = True
    if skip:
        continue
    parts = line.split()
    if len(parts) > 1:
        degauss_value = float(parts[-1])
        N_Ef_Ry_value = float(parts[11])
        N_Ef_eV_value = 2 * N_Ef_Ry_value / Ry2ev
        degauss.append(degauss_value)
        N_Ef_Ry.append(N_Ef_Ry_value)
        N_Ef_eV.append(N_Ef_eV_value)

# 保存数据为 CSV 文件，以 states/eV/Unit Cell (即：states/eV/f.u.) 为单位
df_eV = pd.DataFrame({'degauss': degauss, 'N(Ef) (states/eV/f.u.)': N_Ef_eV})
df_eV.to_csv('elephsigma_TDOS_for_eV.csv', index=False)

# 如果需要还可以保存原始的数据，以states/spin/Ry/Unit Cell (即：states/spin/Ry/f.u.) 为单位。
df_Ry = pd.DataFrame({'degauss': degauss, 'N(Ef) (states/spin/Ry/Unit Cell)': N_Ef_Ry})
df_Ry.to_csv('elephsigma_TDOS_for_Ry.csv', index=False)

# 绘制图像
plt.figure(figsize=(10, 6))
plt.plot(degauss, N_Ef_Ry, marker='o', label='states/Ry/Unit Cell')
plt.xlabel('degauss')
plt.ylabel('N(Ef)')
plt.title('N(Ef) vs degauss')
plt.grid(True)
plt.legend()
plt.savefig('lambda_plot.png')
plt.show()