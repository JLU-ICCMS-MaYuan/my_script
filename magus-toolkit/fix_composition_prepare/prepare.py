#!/usr/bin/env python

import os
import sys

def prepare_files(dirs_path, mode=None, formula=None, minat=None, maxat=None):
    if os.path.exists('prepare'):
        os.system(f'cp -rf prepare/* {dirs_path}')
    inputyaml = os.path.join(dirs_path, 'input.yaml')
    with open(inputyaml) as f:
        lines = f.readlines()
    for i, line in enumerate(lines):
        if mode == "Great_power_produces_miracles":
            if 'formula:' in line:
                lines[i] = 'formula: [[1,0,0],[0,1,0],[0,0,1]]\n'
            if 'min_n_atoms:' in line:
                lines[i] = 'min_n_atoms: 3\n'
            if 'max_n_atoms:' in line:
                lines[i] = 'max_n_atoms: 30\n'
            if 'full_ele:' in line:
                lines[i] = 'full_ele: False\n'
        elif mode == "Divide_and_conquer":
            if 'formula:' in line:
                x1, x2, x3, y1, y2, y3, z1, z2, z3 = formula
                lines[i] = f'formula: [[{x1},{x2},{x3}],[{y1},{y2},{y3}],[{z1},{z2},{z3}]]\n'
            if 'min_n_atoms:' in line :
                lines[i] = f'min_n_atoms: {minat}\n'
            if 'max_n_atoms:' in line:
                lines[i] = f'max_n_atoms: {maxat}\n'
            if 'full_ele:' in line:
                lines[i] = 'full_ele: True\n'
    with open(inputyaml, 'w') as f:
        f.writelines(lines)

         

if __name__ == '__main__':
    info = '''Attention to set these parameters in input.yaml by yourself:
symbols
d_ratio
volume_ratio
MainCalculator/preProcessing
MainCalculator/ppLabel
MainCalculator/numCore
MainCalculator/queueName
MLCalculator/preProcessing
MLCalculator/queueName
MLCalculator/numCore
MLCalculator/min_dist

python prepare.py 200   1 4 0   0 1 1   0 0 1
200: pressure
1 4 1   0 1 1   0 0 1 : formula
'''
    print(info)
    atoms_range = [[3,10],[11,15],[16,20],[21,25],[26,30],[31,35],[35,40]]
    pressure = sys.argv[1]
    formula  = sys.argv[2:]
    
    dirs_path = '{}/{}.{}-{}_{}GPa'.format(pressure, 0, 3, 30, pressure)
    mode="Great_power_produces_miracles"
    print(dirs_path)
    if not os.path.exists(dirs_path):
        os.system(f' mkdir -p {dirs_path}')
    prepare_files(dirs_path, mode)
    
    
    mode="Divide_and_conquer"
    for i, [minat, maxat] in enumerate(atoms_range): 
        dirs_path = '{}/{}.{}-{}_{}GPa'.format(pressure, i+1, minat, maxat, pressure)
        print(dirs_path)
        if not os.path.exists(dirs_path):
            os.system(f' mkdir -p {dirs_path}')
        prepare_files(dirs_path, mode, formula, minat, maxat)
