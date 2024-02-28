#!/usr/bin/env python3
import re
import os
import sys

from pathlib import Path
from pprint import pprint

from ase.io import read


def get_dH(outcar_path):
    try:
        dH = os.popen(f"grep -a enthalpy {outcar_path} " + " | tail -n 1 | awk '{print $ 5}'").read().strip('\n')
        if not dH:
            dH = 100000000000000.0
        else:
            dH = float(dH)
    except:
        print("   OUTCAR有问题读不出来焓值")
        dH = 100000000000000.0
    return dH

def get_convexhull_info(poscar_path, outcar_path):
    atoms = read(poscar_path)
    dH = get_dH(outcar_path)
    N = atoms.get_global_number_of_atoms()
    # empirical=True 表示 获得化学式的最简化学式
    formula = atoms.get_chemical_formula(empirical=True)
    dH_peratom = dH/N
    return formula, dH_peratom


if __name__ == "__main__":
    info='''注意：-------------------------------
        这个脚本必须在压强目录的上一层目录中使用，比如：
            1.LaY
                1.La1-Y1  
                    10 30 50 80 100 150 200
                2.La3-Y1  
                    10 30 50 80 100 150 200
                3.La1-Y3
                    10 30 50 80 100 150 200
            2.LaCe
                ....
            3.YCe
                ....
        你必须在与1.La1-Y1, 2.La3-Y1 ,3.La1-Y3同一级的目录中使用
    
    使用方法是：
        makeConvexHullData.py 压强值
    
    该脚本会自动将1.La1-Y1, 2.La3-Y1, 3.La1-Y3中相同压强值的结构的焓值提取出来
    '''
    print(info)
    press = sys.argv[1]
    fail_d = []
    success_d = []
    none_d = []
    current_directory = Path.cwd()
    directories = [item.name for item in current_directory.iterdir() if item.is_dir()]
    convexhull_file = open(str(press)+"GPa"+"convexhull.csv", "w")
    print("{},{},{}".format("Number","formula","enthalpy"), file=convexhull_file)
    for dirs in directories:
        work_path = Path(dirs).joinpath(str(press))
        outcar_path = work_path.joinpath("OUTCAR")
        poscar_path = work_path.joinpath("POSCAR")
        if outcar_path.exists() and poscar_path.exists():
            idx = work_path.parent.name.split(".")[0]
            # 抑制错误消息
            res = os.popen(f'grep -s "reached required accuracy - stopping structural energy minimisation" {outcar_path}').read() # 如果没有找到指定内容不输出错误结果。
            if res:
                formula, dH_peratom = get_convexhull_info(poscar_path , outcar_path)
                print("{},{},{}".format(idx, formula, dH_peratom), file=convexhull_file)
            else:
                fail_d.append(work_path.__str__())
        else:
            none_d.append(work_path.__str__())
