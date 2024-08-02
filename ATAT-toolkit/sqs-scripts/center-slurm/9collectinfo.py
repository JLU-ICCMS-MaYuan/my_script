#!/usr/bin/env python3
import os
import re
import sys

from pathlib import Path

import pandas as pd

from pymatgen.core.structure import Structure

def get_enthalpy(outcar_path):
    enthalpy_list = []
    with open(outcar_path, "r") as outcar:
        for line in outcar.readlines():
            if "enthalpy" in line:
                enthalpy = re.search(r"\-*\d+\.\d+", line).group()
                enthalpy_list.append(enthalpy)
    return enthalpy_list[-1]

def get_volume(outcar_path):
    v = float(os.popen(f"""grep volume {outcar_path} | tail -n 1 """ + """| awk -F':' '{print $2}'""").read().strip())
    return v

if __name__ == "__main__":
    
    press = input("输入要提取能量的压强\n")
    info = []
    for drt in Path.cwd().glob(r"POSCAR-*"):
        print(drt)
    #    outcar_path = drt.joinpath(str(press),"OUTCAR")
    #    poscar_path = drt.joinpath(str(press),"POSCAR")
        outcar_path = drt.joinpath("OUTCAR")
        poscar_path = drt.joinpath("POSCAR")
        try:
            struct = Structure.from_file(poscar_path)  
            total_atoms = struct.composition.num_atoms 
            enthalpy = get_enthalpy(outcar_path)
            enthalpy_per_atoms = float(enthalpy) / total_atoms
            volume = get_volume(outcar_path)
            info.append([drt.name, enthalpy_per_atoms, volume])
        except:
            print(f"{drt.name} OUTCAR 或者 POSCAR 有问题")

    data = pd.DataFrame(info, columns=["name", "enthalpy(eV/atom)", "volume(A^3)"])
    data.to_csv(str(press)+"GPa-"+"dH.csv", index=None)

