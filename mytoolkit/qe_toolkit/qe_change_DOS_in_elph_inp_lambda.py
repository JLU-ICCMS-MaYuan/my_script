#!/usr/bin/env python3
import sys

irr_num = sys.argv[1]

with open('elph.inp_lambda.1', 'r') as e1:
    e1_lines = e1.readlines()

e1DOSs_EFs = []
for e1line in e1_lines:
    if 'DOS' in e1line and 'Ef=' in e1line:
        e1DOSs_EFs.append(e1line)

if len(e1DOSs_EFs) != 10:
    print(f'{len(e1DOSs_EFs)} 个 DOS值')
    sys.exit(0)


for i in range(2,int(irr_num)+1):
    with open(f'elph.inp_lambda.{i}', 'r') as ei:
        ei_lines = ei.readlines()
    j = 0
    for idx, eiline in enumerate(ei_lines):
        if 'DOS' in eiline and 'Ef=' in eiline:
            ei_lines[idx] = e1DOSs_EFs[j]
            j += 1
    with open(f'elph.inp_lambda.{i}', 'w') as ei:
        ei.writelines(ei_lines)
