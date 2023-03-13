#!/usr/bin/env python3
# code by hxl and her shimei 2023.3.11

import sys
enthalpy_A, enthalpy_B, enthalpy_C = sys.argv[1:4]

print(f"已输入的赝三元相图的端点值 {enthalpy_A} {enthalpy_B} {enthalpy_C}")

with open('Individuals')as f:
	data = [i for i in f.readlines()[2:]]
data = [i for i in data if not 'N/A' in i]
data = [i for i in data if (not 'e-' in i) and (not 'e+' in i)]
data = [i.split() for i in data]
#data = [i for i in data if len(i) in [19, 20]]
#        0            8            46           238.668                 156.773  
data = [[float(i[4]), float(i[5]), float(i[6]), float(i[8])] for i in data]
#[H_MgO, H_SiO2, H_H2O] = [12.9305, 10.996, 3.012]
# [H_Y, H_Sr, H_H] = [20.8371, 27.7716, 1.0010]

result = []
for i in data:
	[n_A, n_B, n_C] = [i[0], i[1], i[2]/2]
	n = n_A + n_B + n_C
	[ratio_A, ratio_B, ratio_C] = [n_A/n, n_B/n, n_C/n]
	
	H = i[3]/n - ratio_A*enthalpy_A - ratio_B*enthalpy_B - ratio_C*enthalpy_C
	result.append([float(format(ratio_A, '.4f')), float(format(ratio_B, '.4f')), float(format(ratio_C, '.4f')), float(format(H, '.4f'))])

lines = ['    x    y    z    H\n']
for i in result:
	i = ['%.4f'%j for j in i]
	lines.append(' '.join(i) + '\n')
with open('out1.dat', 'w')as f:
	f.writelines(lines)

result_d = {}
for i in result:
	[x, y, z] = i[:3]
	if not (x in result_d):
		result_d[x] = {}
	if not (y in result_d[x]):
		result_d[x][y] = {}
	if not (z in result_d[x][y]):
		result_d[x][y][z] = i[3]
	if i[3] < result_d[x][y][z]:
		result_d[x][y][z] = i[3]

lines = ['    x    y    z    H\n']
for x in result_d:
	for y in result_d[x]:
		for z in result_d[x][y]:
			i = [x, y, z, result_d[x][y][z]]
			i = ['%.4f'%j for j in i]
			lines.append(' '.join(i) + '\n')
with open('out2.dat', 'w')as f:
	f.writelines(lines)

