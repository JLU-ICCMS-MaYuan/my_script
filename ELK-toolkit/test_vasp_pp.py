#!/usr/bin/env python3
import os
import sys
import shutil

import numpy as np

def write_incar(pstress, work_path):
    incar_path = os.path.join(work_path, "INCAR")
    with open(incar_path, 'w') as incar:
        incar.write("ISTART   = 0    \n")   
        incar.write("ICHARG   = 2    \n")
        incar.write("ISYM     = 2    \n") 
        incar.write("ENCUT    = 600  \n")
        incar.write("PREC     = A    \n") 
        incar.write("NCORE    = 4    \n")         
        incar.write("KSPACING = 0.188\n")            
        incar.write("ISMEAR   = 0    \n")   
        incar.write("SIGMA    = 0.2  \n")   
        incar.write("NELM     = 200  \n")   
        incar.write("NELMIN   = 6    \n")   
        incar.write("EDIFF    = 1e-6 \n")
        incar.write("EDIFFG   = -0.1 \n")
        incar.write("NSW      = 200  \n")   
        incar.write("IBRION   = 2    \n")   
        incar.write("ISIF     = 3    \n")
        incar.write("POTIM    = 0.05  \n")
        incar.write("LWAVE    = .FALSE.\n")                 
        incar.write("LCHARG   = .FALSE.\n")   
        incar.write("PSTRESS  = {}   \n".format(pstress*10))

def get_totalenergy(outcar_path):
    eV2hartree = 0.036749322175655
    hartree2eV = 27.211386245988
    try:
        dE = float(os.popen(f'grep -s "free  energy   TOTEN" {outcar_path}' + " | tail -n 1 | awk '{print $ 5}'").read().strip('\n'))
        dE = dE * eV2hartree 
    except:
        dE = 100000000000
    return dE

def get_volume(outcar_path):
    try:
        V = float(os.popen(f'grep -s "volume of cell" {outcar_path}' + " | tail -n 1 | awk '{print $ 5}'").read().strip('\n'))
        V = V / np.power(0.529, 3)
    except:
        V = 100000000000
    return V


if __name__ == "__main__":
    
    print("Note: --------------------")
    print("    你需要在当前目录下准备好: POSCAR, POTCAR, INCAR, submit.sh")
    print("    测试的PSTRESS值分别是: -5 0.0001 50 100 150 200 250 300 350 400")
    print("    该脚本不提供自动提任务的命令: 你可以用以下命令提供命令:")
    print("        for i in p-5 p+0.00001 p+50 p+100 p+150 p+200 p+250 p+300 p+350 p+400; do cd $i; qsub submit.sh; cd ..; done")
    print("        for i in p-5 p+0.00001 p+50 p+100 p+150 p+200 p+250 p+300 p+350 p+400; do cd $i; sbatch submit.sh; cd ..; done")
    print("\nNote: --------------------")
    print("    创建测试VASP的SIGMA输入文件目录以及准备vasp的输入文件")

    press_s = [-5, 0.00001, 50, 100, 150, 200, 250, 300, 350, 400]
    potcar_path = os.path.abspath("POTCAR")
    poscar_path = os.path.abspath("POSCAR")
    submit_path = os.path.abspath("submit.sh")

    for pstress in press_s:
        if type(pstress) == int :
            test_path = os.path.abspath('p{:+d}'.format(pstress))
        elif type(pstress) == float:
            test_path = os.path.abspath('p{:+.5f}'.format(pstress))
            
        if not os.path.exists(test_path):
            os.mkdir(test_path)
        os.system(f"cp -f {potcar_path} {test_path}")
        os.system(f"cp -f {poscar_path} {test_path}")
        os.system(f"cp -f {submit_path} {test_path}")
        write_incar(pstress, test_path)

    V_pstress = []
    v_energy = []
    for pstress in press_s:
        if type(pstress) == int :
            outcar_path = os.path.join('p{:+d}'.format(pstress), "OUTCAR")
        elif type(pstress) == float:
            outcar_path = os.path.join('p{:+.5f}'.format(pstress), "OUTCAR")
        E = get_totalenergy(outcar_path)
        V = get_volume(outcar_path)
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
        print("All OUTCARs are OK, V_pstress.csv has been wroten in current position")
    else:
        print("If all OUTCARs are OK, V_pstress.csv will be wroten in current position")

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
        print("All OUTCARs are OK, V_pstress.csv has been wroten in current position")
    else:
        print("If all OUTCARs are OK, V_pstress.csv will be wroten in current position")

    print("Note: --------------------")
    print("    现在获得了eos.in文件, 直接./eos就可以获得PVPAI.OUT文件了")
    print("    HPPAI.OUT: (H-P)的关系")
    print("    PARAM.OUT: 体弹模量")
