#!/usr/bin/env python3
import os
import sys
import shutil

import numpy as np




def get_totalenergy(scf_path):
    try:
        dE = float(os.popen(f'grep -s "TOTAL ENERGY IN Ry" {scf_path}' + " | tail -n 1 | awk '{print $ 9}'").read().strip('\n'))
        dE = dE/2
    except:
        dE = 100000000000
    return dE

def get_volume(scf_path):
    try:
        V = float(os.popen(f'grep -s "UNIT CELL VOLUME" {scf_path}' + " | tail -n 1 | awk '{print $ 7}'").read().strip('\n'))
    except:
        V = 100000000000
    return V


if __name__ == "__main__":
    pstress_s = [0.00001, 50, 100, 150, 200, 250, 300, 350, 400]

    for pstress in pstress_s:
        if type(pstress) == int :
            src_file = os.path.join('p{:+d}'.format(pstress)+'.cif')
            test_path = os.path.abspath('p{:+d}'.format(pstress))
            dst_file = os.path.join(test_path, 'p{:+d}'.format(pstress)+'.cif')
        elif type(pstress) == float:
            src_file = os.path.join('p{:+.5f}'.format(pstress)+'.cif')
            test_path = os.path.abspath('p{:+.5f}'.format(pstress))
            dst_file = os.path.join(test_path, 'p{:+.5f}'.format(pstress)+'.cif')
        if not os.path.exists(test_path):
            os.mkdir(test_path)
        shutil.copy(src_file, dst_file)

    print("Note: --------------------")
    print("    你需要在最高压强400GPa下计算wien2k的scf, 然后获得RMT值。")
    print("    然后修改其余压强下的case.struct中的RMT值为400GPa下的RMT值。")


    V_pstress = []
    v_energy = []
    for pstress in pstress_s:
        if type(pstress) == int :
            scf_path = os.path.join('p{:+d}'.format(pstress), 'p{:+d}'.format(pstress)+'.scf')
        elif type(pstress) == float:
            scf_path = os.path.join('p{:+.5f}'.format(pstress), 'p{:+.5f}'.format(pstress)+'.scf')
        E = get_totalenergy(scf_path)
        V = get_volume(scf_path)
        V_pstress.append([V, pstress])
        v_energy.append([V, E])

    V_pstress = np.array(V_pstress)
    if len(V_pstress) == 9:
        print("{:<14},{:<14}".format("V(A^3)", "pstress(GPa)"))
        with open("V_pstress.csv", 'w') as f:
            f.write("{:<14},{:<14}\n".format("V(A^3)", "pstress(GPa)"))
            for V, pstress in V_pstress:
                f.write("{:<14.8f},{:<14.8f}\n".format(V, pstress))
                print("{:<14.8f},{:<14.8f}".format(V, pstress))
        print("All case.scf are OK, V_pstress.csv has been wroten in current position")
    else:
        print("If all case.scf are OK, V_pstress.csv will be wroten in current position")

    v_energy = np.array(v_energy)
    if len(v_energy) == 9:
        print("{:<14},{:<14}".format("V(bohr^3)", "E(hartree)", ))
        with open("eos.in", 'w') as f:
            name             = input("晶体结构名称(建议:crystal)\n");             f.write("{}\n".format(name))
            numatoms         = input("对应原子总数目(建议:胞内原子数)\n");         f.write("{}\n".format(numatoms))
            bm_type          = input("选择的是第几种物态方程(共7种, 建议:3)\n");   f.write("{}\n".format(bm_type))
            insert_points    = input("拟合时每两对E-V之间的插点数目(建议:1000)\n");f.write("{}  {}  {}\n".format(v_energy[0][0], v_energy[-1][0], insert_points))
            num_press_points = input("共有多少个压力点(自己数)\n");               f.write("{}\n".format(num_press_points))
            for V, E in v_energy:
                f.write("{:<14.8f}  {:<14.8f}\n".format(V, E))
                print("{:<14.8f}  {:<14.8f}".format(V, E))
        print("All case.scf are OK, V_pstress.csv has been wroten in current position")
    else:
        print("If all case.scf are OK, V_pstress.csv will be wroten in current position")


