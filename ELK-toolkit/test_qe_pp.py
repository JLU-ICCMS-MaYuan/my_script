#!/usr/bin/env python3

import os
import sys

import numpy as np

def get_cell(poscar_path):
    sed_order = "sed -n '3,5p' " + f" {poscar_path} "
    content = os.popen(sed_order).readlines()
    cell_parameters = [cell.strip("\n") for cell in content]
    return cell_parameters

def get_coords(poscar_path):
    sed_order1 = "sed -n '6p' " + f" {poscar_path} "
    elements   = os.popen(sed_order1).read().strip('\n').split()
    sed_order2 = "sed -n '7p' " + f" {poscar_path} "
    numatoms   = os.popen(sed_order2).read().strip('\n').split()
    # 检查元素的个数与原子种类的个数是否相同
    assert len(elements) == len(numatoms)
    # 检查该POSCAR是否是分数坐标
    sed_order3 = "sed -n '8p' " + f" {poscar_path} "
    coordstype = os.popen(sed_order3).read().strip('\n')
    assert coordstype[0] == 'D' or coordstype[0] == 'd'
    # 将元素按照原子个数重复写进列表elementlist里
    sed_order4 = "sed -n '9,$p' " + f" {poscar_path} "
    coordnates = os.popen(sed_order4).readlines()
    element_names= []
    for na, ele in zip(numatoms, elements): # extend() 函数用于在列表末尾一次性追加另一个序列中的多个值
        element_names.extend([ele]*int(na))

    # 建立元素列表elementlist与分数坐标的一一对应关系
    fractional_sites = list(zip(element_names, coordnates))
    # 按照elements的元素顺序去写原子坐标，因为赝势的顺序就是elements的顺序
    # 这样就可以保证赝势的顺序pp_order和原子坐标的顺序一致！！！！
    pp_order = [spe for spe in elements]
    fractional_sites = sorted(fractional_sites, key=lambda item: pp_order.index(item[0]))
    fractional_sites = [
        '{:<4}   {}'.format(ele, site.strip("\n")) for ele, site in fractional_sites]
    return fractional_sites

def change_info(old_scf_in_path, new_scf_in_path, press, poscar_path=None):

    with open(old_scf_in_path, "r") as f:
        lines = f.readlines()

    for idx, line in enumerate(lines):
        if "press = " in line:
            lines[idx] = f" press = {press*10}\n"
        if "pseudo_dir" in line:
            lines[idx] = " pseudo_dir='../pp'\n"
        if "CELL_PARAMETERS" in line:
            cell_idx = idx
        if "ATOMIC_POSITIONS" in line:
            coords_idx = idx

    if poscar_path:
        cell = get_cell(poscar_path)
        coords = get_coords(poscar_path)
        for i, cell_p in enumerate(cell):
            lines[cell_idx+i+1] = cell_p+'\n'
        for i, coord in enumerate(coords):
            lines[coords_idx+i+1] = coord+'\n'

    with open(new_scf_in_path, "w") as f:
        f.writelines(lines)

def get_total_energy(scf_out_path):
    """计算并返回结构能量 单位的Ry"""
    try:
        Et = float(os.popen(f""" grep -s "!    total energy" {scf_out_path}""" + """ | tail -n 1  | awk '{print $5}' """).read().strip('\n'))
        Et = Et/2 # 注意单位是Ry
        return Et
    except:
        Et = 1000000000000.0
        return Et

def get_volume(scf_out_path):
    angstrom2bohr = 0.14803588900000003
    bohr2angstrom = 6.75511868611806558
    try:
        # V = float(os.popen(f'grep -s "new unit-cell volume" {scf_out_path}' + " | tail -n 1 | awk '{print $ 5}'").read().strip('\n'))
        V = float(os.popen(f'grep -s "unit-cell volume" {scf_out_path}' + " | tail -n 1 | awk '{print $ 4}'").read().strip('\n'))
        # V = V * (0.52917721092**3)
    except:
        V = 100000000000
    return V

if __name__ == "__main__":

    print("Note: --------------------")
    print("    你需要在当前目录下准备好: scf.in, pp, submit.sh(pp目录中放好赝势)")
    print("    测试的press值分别是: -5 0.0001 50 100 150 200 250 300 350 400")
    print("    该脚本不提供自动提任务的命令: 你可以用以下命令提供命令:")
    print("        for i in p-5 p+0.00001 p+50 p+100 p+150 p+200 p+250 p+300 p+350 p+400; do cd $i; qsub submit.sh; cd ..; done")
    print("        for i in p-5 p+0.00001 p+50 p+100 p+150 p+200 p+250 p+300 p+350 p+400; do cd $i; sbatch submit.sh; cd ..; done")

    print("Note: --------------------")
    print("    创建测试qe的赝势输入文件目录以及准备qe的输入文件")
    press_s = [-5, 0.00001, 50, 100, 150, 200, 250, 300, 350, 400]
    old_scf_in_path = os.path.abspath("scf.in")
    for press in press_s:
        poscar_path = None
        if type(press) == int :
            test_path = os.path.abspath('p{:+d}'.format(press))
            new_scf_in_path = os.path.join(test_path, "scf.in")
            if sys.argv[1]:
                poscar_path = os.path.join(sys.argv[1], 'p{:+d}'.format(press),"CONTCAR")
        elif type(press) == float:
            test_path = os.path.abspath('p{:+.5f}'.format(press))
            new_scf_in_path = os.path.join('p{:+.5f}'.format(press), "scf.in")
            if sys.argv[1]:
                poscar_path = os.path.join(sys.argv[1], 'p{:+.5f}'.format(press),"CONTCAR")
        if not os.path.exists(test_path):
            os.mkdir(test_path)
        change_info(old_scf_in_path, new_scf_in_path, press, poscar_path)
        submit_path = os.path.abspath("submit.sh")
        os.system(f"cp -f {submit_path} {test_path}")




    V_pstress = []
    v_energy = []
    for press in press_s:
        if type(press) == int :
            test_path = os.path.abspath('p{:+d}'.format(press))
            scf_out_path = os.path.join(test_path, "scf.out")
        elif type(press) == float:
            test_path = os.path.abspath('p{:+.5f}'.format(press))
            scf_out_path = os.path.join('p{:+.5f}'.format(press), "scf.out")
        E = get_total_energy(scf_out_path)
        V = get_volume(scf_out_path)
        V_pstress.append([V, press])
        v_energy.append([V, E])

    V_pstress = np.array(V_pstress)
    if len(V_pstress) == 9:
        print("{:<14},{:<14}".format("V(bohr^3)", "pstress(GPa)"))
        with open("V_pstress.csv", 'w') as f:
            f.write("{:<14},{:<14}\n".format("V(A^3)", "pstress(GPa)"))
            for V, pstress in V_pstress:
                f.write("{:<14.8f},{:<14.8f}\n".format(V, pstress))
                print("{:<14.8f},{:<14.8f}".format(V, pstress))
        print("All scf.out are OK, V_pstress.csv has been wroten in current position")
    else:
        print("If all scf.out are OK, V_pstress.csv will be wroten in current position")

    v_energy = np.array(v_energy)
    if len(v_energy) == 10:
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
        print("All scf.out are OK, V_pstress.csv has been wroten in current position")
    else:
        print("If all scf.out are OK, V_pstress.csv will be wroten in current position")

    print("Note: --------------------")
    print("    现在获得了eos.in文件, 直接./eos就可以获得PVPAI.OUT文件了")
    print("    HPPAI.OUT: (H-P)的关系")
    print("    PARAM.OUT: 体弹模量")
