#!/usr/bin/env python3
import re
import os
import sys

from pathlib import Path
from pprint import pprint

from ase.io import read


def get_dH(outcar_path):
    try:
        dH = os.popen(f"grep enthalpy {outcar_path} " + " | tail -n 1 | awk '{print $ 5}'").read().strip('\n')
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

    press = sys.argv[1]
    fail_d = []
    success_d = []
    none_d = []
    current_directory = Path.cwd()
    directories = [item.name for item in current_directory.iterdir() if item.is_dir()]
    convexhull_file = open("convexhull.csv", "w")
    print("{},{},{}".format("Number","formula","enthalpy"), file=convexhull_file)
    for dirs in directories:
        work_path = Path(dirs).joinpath(str(press))
        ourcar_path = work_path.joinpath("OUTCAR")
        poscar_path = work_path.joinpath("POSCAR")
        if ourcar_path.exists() and poscar_path.exists():
            idx = work_path.parent.name.split("-")[0]
            # 抑制错误消息
            res = os.popen(f'grep -s "reached required accuracy - stopping structural energy minimisation" {ourcar_path}').read() # 如果没有找到指定内容不输出错误结果。
            if res:
                formula, dH_peratom = get_convexhull_info(poscar_path , ourcar_path)
                print("{},{},{}".format(idx, formula, dH_peratom), file=convexhull_file)
            else:
                fail_d.append(work_path.__str__())
        else:
            none_d.append(work_path.__str__())
