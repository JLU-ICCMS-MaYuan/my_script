#!/usr/bin/env python3
import os
import sys
import shutil

import pandas as pd
import numpy as np

from ase.calculators.calculator import kpts2mp
from ase.io import read

import matplotlib.pyplot as plt


def write_elk_in(rgkmax, elk_in, work_path):
   
    with open(elk_in, 'r+') as f:
        lines = f.readlines()
        # 寻找ngridk关键字并修改下一行的值
        for j in range(len(lines)):
            if 'rgkmax' in lines[j].lower():
                lines[j + 1] = f'{rgkmax}\n'

    new_elk_in = os.path.join(work_path, "elk.in")
    with open(new_elk_in, "w") as f:
        f.writelines(lines)

def get_info_per_atom(atom, totenergy_out_path):

    try:
        natoms = atom.get_global_number_of_atoms()
    except:
        natoms = None

    try:
        energy = float(os.popen(f'tail -n 1 {totenergy_out_path}').read().strip('\n'))*27.211386245988
    except:
        energy = None

    # try:
    #     force = np.array(list(map(float, os.popen("grep -A%d 'TOTAL-FORCE' %s |tail -n %d|awk '{print $4,$5,$6}'"%(natoms+1,outcar_path,natoms)).read().split())))
    # except:
    #     force = None

    # try:
    #     virial = np.array(list(map(float, os.popen("grep -A20 '\-STRESS' %s |grep Total|awk '{print $2,$3,$4,$5,$6,$7}'"%(outcar_path)).read().split())))
    # except:
    #     virial = None
  
    if (natoms is not None) and (energy is not None): #and (force is not None) and (virial is not None):
        dE_per_atom = np.array(energy)/natoms
        # virial_per_atom = np.array(virial)/natoms
    else:
        dE_per_atom = 10000000000000000
        # force = 10000000000000000 
        # virial_per_atom = 10000000000000000
    return dE_per_atom #, force, virial_per_atom



if __name__ == "__main__":
    
    # jobsystem = "chaoyin_pbs"
    jobsystem = "tangB_slurm"
    # jobsystem = "coshare_slurm"
    print("Note: --------------------")
    print("    你需要在当前目录下准备好: POSCAR, elk.in, submit.sh")
    print("Note: --------------------")
    start, stop, step = sys.argv[1:]
    print("    python test_vasp_rgkmax.py  1.0 10.0  1.0")
    print("    python test_vasp_rgkmax.py 10.0  1.0 -1.0")
    print("    以上两种都是产生1到10, 间隔为1的rgkmax")
    rgkmaxs = list(np.round(np.arange(eval(start), eval(stop), eval(step)), decimals=2))
#    except:
#        print("你设置的有问题，采用一下的k点密度")
#        kresolutions = [0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.19, 0.18, 0.17, 0.16, 0.15, 0.14, 0.13, 0.12, 0.1, 0.05]
    print("Note: --------------------")
    print("    该脚本不提供自动提任务的命令: 你可以用以下命令提供命令:")
    print("        for i in {}; do cd $i; qsub submit.sh;   cd ..; done".format(' '.join(list(map(str, rgkmaxs)))))
    print("        for i in {}; do cd $i; sbatch submit.sh; cd ..; done".format(' '.join(list(map(str, rgkmaxs)))))

    poscar_path = os.path.abspath("POSCAR")
    elk_in_path  = os.path.abspath("elk.in")
    submit_path  = os.path.abspath("submit.sh")
    
    atom = read(poscar_path)
    
    for rgkmax in rgkmaxs:
        test_path = os.path.abspath(str(rgkmax))
        if not os.path.exists(test_path):
            os.mkdir(test_path)
        write_elk_in(rgkmax, elk_in_path, test_path)
        os.system(f"cp -f {submit_path} {test_path}")

    print("    尝试获得每原子的焓值 dE(meV/atom)")
    rgkmax_dE_dF_dV = []
    for rgkmax in rgkmaxs:
        totenergy_out_path = os.path.join(str(rgkmax), "TOTENERGY.OUT")
        # dE_per_atom, force, virial_per_atom = get_info_per_atom(atom, totenergy_out_path)
        # rgkmax_dE_dF_dV.append([rgkmax, dE_per_atom, force, virial_per_atom])
        dE_per_atom = get_info_per_atom(atom, totenergy_out_path)
        rgkmax_dE_dF_dV.append([rgkmax, dE_per_atom])
    rgkmax_dE_dF_dV = sorted(rgkmax_dE_dF_dV, key=lambda x: x[0],  reverse=False)

    # 方式一：通过相邻两个rgkmax的能量的差判断收敛性
    # Ediff = np.diff(rgkmax_dE, axis=0)[:,-1]
    # Ediff = np.insert(Ediff, 0, [0.0])
    # Ediff = Ediff[:, np.newaxis]
    # rgkmax_dE_Hdiff = np.hstack((rgkmax_dE, Ediff))

    # 方式二：所有kpoints的能量减去最密rgkmax的能量判断收敛性

    with open("rgkmax_dE.csv", 'w') as f:
        f.write("{:<12},{:<14},{:<14},{:<14}\n".format("rgkmax", "Ediff(meV/atom)", "Fdiff(meV/A)", "Vdiff(meV/atom)"))
        print("{:<12},{:<14},{:<14},{:<14}".format("rgkmax", "Ediff(meV/atom)", "Fdiff(meV/A)", "Vdiff(meV/atom)"))
        for idx, kefv in enumerate(rgkmax_dE_dF_dV):
            kmeshes = kefv[0] 
            delta_E = kefv[1] - rgkmax_dE_dF_dV[-1][1]
            # delta_F = kefv[2] - rgkmax_dE_dF_dV[-1][2]
            # delta_V = kefv[3] - rgkmax_dE_dF_dV[-1][3]
            # f.write("{:<12.3f},{:<14.8f},{:<14.8f},{:<14.8f}\n".format(kmeshes, delta_E*1000, np.nanmax(delta_F)*1000,np.nanmax(delta_V)*1000))
            # print("{:<12.3f},{:<14.8f},{:<14.8f},{:<14.8f}".format(kmeshes, delta_E*1000, np.nanmax(delta_F)*1000, np.nanmax(delta_V)*1000))
            f.write("{:<12.3f},{:<14.8f}\n".format(kmeshes, delta_E*1000))
            print("{:<12.3f},{:<14.8f}".format(kmeshes, delta_E*1000))
        else:
            print("If all elk are OK, rgkmax_dE.csv will be wroten in current position")

    # 绘制图表
    df = pd.read_csv("rgkmax_dE.csv")
    plt.figure(figsize=(10, 6))  # 设置图形大小
    plt.plot(df[rgkmax], df['Ediff(meV/atom)'], marker='o', linestyle='-', color='r', label='Ediff (meV/atom)')
    plt.title('rgkmax vs Total Energy')  # 设置图表标题
    plt.xlabel('rgkmax')  # 设置x轴标签
    plt.ylabel('Total Energy (Hartree)')  # 设置y轴标签
    plt.grid(True)  # 显示网格
    plt.savefig('elk_rgkmax_conv.png')  # 保存图形
    plt.show()  # 显示图形
