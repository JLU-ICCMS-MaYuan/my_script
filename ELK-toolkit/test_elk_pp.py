#!/usr/bin/env python
import os
import argparse
import numpy as np

from ase.io import read
from ase.units import Bohr, Hartree
from ase.io.elk import read_elk, write_elk_in
from ase.calculators.calculator import kpts2mp

def get_ngridk_vkloff(atoms, kspacing):
    ngridk = kpts2mp(atoms, kspacing)
    vkloff = []  # is this below correct?
    for nk in ngridk:
        if nk % 2 == 0:  # shift kpoint away from gamma point
            vkloff.append(0.5)
        else:
            vkloff.append(0)
    return ngridk, vkloff

def write_elk_in(input_file, output_file, species_path, kspacing):
    # 读取结构
    atoms = read(input_file)
    task = 0
    latvopt = 1
    epspot = 1.e-7
    epsstress = 1.e-3
    mixtype = 3 # type of mixing required for the potential   3 Broyden mixing
    xctype = 20 # GGA, Perdew-Burke-Ernzerhof, Phys. Rev. Lett. 77, 3865 (1996)
    isgkmax = -2 # species for which the muffin-tin radius will be used for calculating gkmax
    rgkmax = 8
    gmaxvr = 18
    ngridk, vkloff = get_ngridk_vkloff(atoms, kspacing)


    with open(output_file, 'w') as fd:
        fd.write('tasks\n  {}\n\n'.format(task))
        fd.write('latvopt\n  {}\n\n'.format(latvopt))
        fd.write('epspot\n  {}\n\n'.format(epspot))
        fd.write('epsstress\n  {}\n\n'.format(epsstress))
        fd.write('mixtype\n  {}\n\n'.format(mixtype))
        fd.write('xctype\n  {}\n\n'.format(xctype))
        fd.write('isgkmax\n  {}\n\n'.format(isgkmax))
        fd.write('rgkmax\n  {}\n\n'.format(rgkmax))
        fd.write('gmaxvr\n  {}\n\n'.format(gmaxvr))
        fd.write('ngridk\n  {}  {}  {}\n\n'.format(ngridk[0], ngridk[1], ngridk[2]))

        fd.write("sppath\n  '{}/'\n\n".format(species_path))
        # cell
        fd.write('avec\n')
        for vec in atoms.cell:
            fd.write('%.14f %.14f %.14f\n' % tuple(vec / Bohr))
        fd.write('\n')
        # atoms
        species = {}
        symbols = []
        for a, (symbol, m) in enumerate(
            zip(atoms.get_chemical_symbols(),
                atoms.get_initial_magnetic_moments())):
            if symbol in species:
                species[symbol].append((a, m))
            else:
                species[symbol] = [(a, m)]
                symbols.append(symbol)
        fd.write('atoms\n%d\n' % len(species))
        # scaled = atoms.get_scaled_positions(wrap=False)
        scaled = np.linalg.solve(atoms.cell.T, atoms.positions.T).T
        for symbol in symbols:
            fd.write(f"'{symbol}.in' : spfname\n")
            fd.write('%d\n' % len(species[symbol]))
            for a, m in species[symbol]:
                fd.write('%.14f %.14f %.14f 0.0 0.0 %.14f\n' %
                        (tuple(scaled[a]) + (m,)))
                
def get_total_energy(info_path):
    """计算并返回结构能量 单位的hartree"""
    try:
        Et = float(os.popen(f""" grep -s "total energy" {info_path}""" + """ | tail -n 1  | awk '{print $4}' """).read().strip('\n'))
        Et = Et*27.211386230 # 注意单位是Hartree, 需要转化为eV
        return Et
    except:
        Et = 1000000000000.0
        return Et

def get_volume(info_path):
    try:
        V = float(os.popen(f""" grep -s "Unit cell volume" {info_path}""" + """ | tail -n 1  | awk '{print $5}' """).read().strip('\n'))
        # V = V * (0.52917721092**3) # 注意单位是Hartree, 需要转化为eV
        return V
    except:
        V = 1000000000000.0
        return V

if __name__ == "__main__":
    info = """Note: --------------------
    You can use it by:
    python test_elk_pp.py -i ../1.vasp_Ce_GW/ -s ~/soft/elk-10.0.15/species -k 0.18

Note: --------------------
    All units are atomic in ELK (Hartree, Bohr, etc.)

    你需要在当前目录下准备好: submit.sh
    测试的press值分别是: -5 0.0001 50 100 150 200 250 300 350 400
    该脚本不提供自动提任务的命令: 你可以用以下命令提供命令:
    for i in p-5 p+0.00001 p+50 p+100 p+150 p+200 p+250 p+300 p+350 p+400; do cd $i; qsub submit.sh; cd ..; done
    for i in p-5 p+0.00001 p+50 p+100 p+150 p+200 p+250 p+300 p+350 p+400; do cd $i; sbatch submit.sh; cd ..; done
Note: --------------------
    创建测试elk的eos输入文件目录以及准备elk的输入文件"""
    print(info)
    parser = argparse.ArgumentParser(description='Generate ELK input file for structure optimization.')
    parser.add_argument('-i', '--input', required=True, help='Input structure directory')
    parser.add_argument('-s', '--sppath', required=True, help='Species path directory')
    parser.add_argument('-k', '--kspacing', type=float, required=True, help='K-point spacing for the grid')

    args = parser.parse_args()
    

    press_s = [-5, 0.00001, 50, 100, 150, 200, 250, 300, 350, 400]
    for press in press_s:
        if type(press) == int:
            test_path = os.path.abspath('p{:+d}'.format(press))
            new_elk_in = os.path.join(test_path, "elk.in")
            poscar_path = os.path.join(args.input, 'p{:+d}'.format(press), "CONTCAR")
        elif type(press) == float:
            test_path = os.path.abspath('p{:+.5f}'.format(press))
            new_elk_in = os.path.join(test_path, "elk.in")
            poscar_path = os.path.join(args.input, 'p{:+.5f}'.format(press),"CONTCAR")

        if not os.path.exists(test_path):
            os.mkdir(test_path)

        write_elk_in(poscar_path, new_elk_in, args.sppath, args.kspacing)
        submit_path = os.path.abspath("submit.sh")
        os.system(f"cp -f {submit_path} {test_path}")

    V_pstress = []
    v_energy = []
    for press in press_s:
        if type(press) == int :
            test_path = os.path.abspath('p{:+d}'.format(press))
            info_out_path = os.path.join(test_path, "INFO.OUT")
        elif type(press) == float:
            test_path = os.path.abspath('p{:+.5f}'.format(press))
            info_out_path = os.path.join('p{:+.5f}'.format(press), "INFO.OUT")
        E = get_total_energy(info_out_path)
        V = get_volume(info_out_path)
        V_pstress.append([V, press])
        v_energy.append([V, E])

    V_pstress = np.array(V_pstress)
    if len(V_pstress) == 9:
        print("{:<14},{:<14}".format("V(A^3)", "pstress(GPa)"))
        with open("V_pstress.csv", 'w') as f:
            f.write("{:<14},{:<14}\n".format("V(A^3)", "pstress(GPa)"))
            for V, pstress in V_pstress:
                f.write("{:<14.8f},{:<14.8f}\n".format(V, pstress))
                print("{:<14.8f},{:<14.8f}".format(V, pstress))
        print("All INFO.OUT are OK, V_pstress.csv has been wroten in current position")
    else:
        print("If all INFO.OUT are OK, V_pstress.csv will be wroten in current position")

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
        print("All INFO.OUT are OK, V_pstress.csv has been wroten in current position")
    else:
        print("If all INFO.OUT are OK, V_pstress.csv will be wroten in current position")

    print("Note: --------------------")
    print("    现在获得了eos.in文件, 直接./eos就可以获得PVPAI.OUT文件了")
    print("    HPPAI.OUT: (H-P)的关系")
    print("    PARAM.OUT: 体弹模量")
