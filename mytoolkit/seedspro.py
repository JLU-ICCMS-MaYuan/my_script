#!/usr/bin/env python3
import os
import re
import sys

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

def read_enthalpy(dir_vr:str, endpoints:list):

    result = []
    for root, dirs, files in os.walk(dir_vr):
        if "OUTCAR" in files and "POSCAR" in files and "CONTCAR" in files:
            print(f"Read in {root}")
            outcar_path = os.path.join(root, "OUTCAR")
            contcar_path= os.path.join(root, "CONTCAR")
            struct = Structure.from_file(contcar_path)

            atoms_amount = struct.composition.num_atoms
            hcontent = float(struct.composition.get_el_amt_dict()['H']) / float(atoms_amount)

            # 计算焓值时，要考虑分子式可能时多倍的。既然取了最简形式的分子式，那么对于的焓值就要除以分子式倍数
            formula = struct.composition.get_integer_formula_and_factor()[0]
            factor  = struct.composition.get_integer_formula_and_factor()[1]
            enthalpy = float(get_enthalpy(outcar_path=outcar_path)) / float(factor)

            number   = "seeds-"+os.path.basename(root)

            series = pd.Series({
                "Number": number,
                "formula": formula,
                "hcontent": hcontent,
                "enthalpy": enthalpy,
            })
            result.append(series)
            

    df = pd.DataFrame(result).sort_values(by='formula')  # 分隔符的用法: \s表示由空格作为分隔符, +表示有多个空格

    print("endpoints单独写一个csv文件----endpoints.csv")
    df_endpoints = df[df["formula"].isin(endpoints)]
    df_endpoints = df_endpoints.sort_values(by="enthalpy").groupby("formula").first().reset_index()
    df_endpoints = df_endpoints.reindex(columns=["Number", "formula", "hcontent", "enthalpy"])
    df_endpoints.to_csv("endpoints.csv", index=False)
    print("Done")

    print("endpoints + 一些之前算的结构单独写一个csv文件----endpoints.csv")
    df_rest = df[~df["formula"].isin(endpoints)]
    df_all = pd.concat([df_endpoints, df_rest], axis=0) # axis=0 表示沿着行的方向进行连接
    df_all.to_csv("seeds.csv", index=False)
    print("Done")


if __name__ == "__main__":
    # opt_dir = sys.argv[1]
    # endpoints = sys.argv[2:]
    print("你必须制定2个输入参数")
    dir_vr = input("结构优化的文件所在的目录路径\n")
    endpoints = input("请指定相图的端点(比如: Lu N2 H2), 只可以使用空格间隔你的端点名称。确保你使用的端点名称在结构文件中存在\n").split()
    read_enthalpy(dir_vr, endpoints)