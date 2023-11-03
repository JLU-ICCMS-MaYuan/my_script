#!/usr/bin/env python3

# '''
# 输出origin绘制convex hull 的数据
# plot_convexhull.py -i ./nnconvexhull.csv -ebh 100 -hand Mg B H
# 保存图片
# plot_convexhull.py -i ./nnconvexhull.csv -ebh 100 -save
# 展示图片
# plot_convexhull.py -i ./nnconvexhull.csv -ebh 100 -show
# 回收所有稳定结构
# plot_convexhull.py -i ./nnconvexhull.csv -ebh 100 -cs -cu
# '''

import sys
import os
import re
import shutil
from pathlib import Path
from argparse import ArgumentParser, RawTextHelpFormatter

import pandas as pd
from pymatgen.analysis.phase_diagram import PDEntry, PhaseDiagram, PDPlotter
from pymatgen.core.composition import Composition
from pymatgen.core.periodic_table import Element


parser = ArgumentParser(formatter_class=RawTextHelpFormatter)
parser.add_argument(
    "-i",
    "--input",
    type=str,
    default="./",
    dest="input_csv_path",
    help="请输入csv文件路径\n"
        "这个csv文件是指: nnconvexhull.csv\n"
        "如果你没有这个文件\n"
        "   请运行`cak3.py --vasp 得到所有的结构搜索到的文件`\n"
        "   然后运行`data_processer.py 得到nnconvexhull.csv 和 得到 nconvexhull.csv`\n"
        "\n"
        "整个命令在执行命令后, 将会在屏幕上输出高于convex hull 0 ~ EnthalpyAboveHullValue meV 的结构 !!!"
)
parser.add_argument(
    "-ed",
    "--endnotes",
    action="store",
    default=None,
    dest="endnotes",
    nargs="+",
    help="如果使用该参数, 需要用户指定某几个单质作为端点值"
)
parser.add_argument(
    "-cs",
    "--collect-stable",
    action="store_true",
    default=False,
    dest="collect_stable",
    help="是否收集稳定结构数据\n"
        "将稳定的结构(落在凸包图)上的结构提取出来放在一个叫`stable_structs`的目录里\n"
)
parser.add_argument(
    "-cu",
    "--collect-unstable",
    action="store_true",
    default=False,
    dest="collect_unstable",
    help="是否收集亚稳结构数据\n"
        "将亚稳的结构(落在凸包图)上的结构提取出来放在一个叫`unstable_structs`的目录里\n"
)
parser.add_argument(
    "-ebh",
    type=float,
    dest="EnthalpyAboveHullValue",
    default=[0, 10.0],
    nargs="+",
    help="高于convex hull xxx emV 的能量的上限\n"
        "在 ebh_lower_limit ~ ebh_higher_limit 这个范围内的亚稳的结构确定出来\n"
        "所以说要输入两个值，第一个值是下陷，第二个值是上限。\n"
)

args = parser.parse_args()

input_csv_path = args.input_csv_path
collect_stable = args.collect_stable
collect_unstable = args.collect_unstable
ebh_lower_limit = args.EnthalpyAboveHullValue[0]
ebh_higher_limit = args.EnthalpyAboveHullValue[1]
endnotes = args.endnotes

# 生成 凸包图对象
# convexhull_data = pd.read_csv(input_csv_path, header=0, sep=',') #  header表示第一行为标题行
print("读入文件中的能量和化学式 (注意：能量必须是化学式的能量，不是每原子的能量) ")
try:
    # 该情况处理的文件：
    # Number  formula     enthalpy
    # 1       Ax1By1Cz1   -1.34343
    # 2       Ax2By2C2z   -2.324324
    convexhull_data = pd.read_table(input_csv_path, header=0, sep=',') #  header表示第一行为标题行
    ini_entries = []
    for idx, row in convexhull_data.iterrows():
        comp = Composition(row['formula'])
        num_at = comp.num_atoms
        enth = row['enthalpy']*num_at
        entry_id = row['Number']
        _entry = PDEntry(comp, enth)
        _entry.entry_id = entry_id
        ini_entries.append(_entry)

except:
    # 该情况处理的文件：
    # Number,formula,enthalpy
    # 1,Ax1By1Cz1,-1.34343
    # 2,Ax2By2C2z,-2.324324
    convexhull_data = pd.read_csv("nnconvexhull.csv", header=0, sep=',')
    ini_entries = []
    for idx, row in convexhull_data.iterrows():
        comp = Composition(row['formula'])
        num_at = comp.num_atoms
        enth = row['enthalpy']*num_at
        entry_id = row['Number']
        _entry = PDEntry(comp, enth)
        _entry.entry_id = entry_id
        ini_entries.append(_entry)

elements_endnotes = [Element(ed) for ed in endnotes]
ini_pd = PhaseDiagram(ini_entries, elements=elements_endnotes)

# 收集所有的稳定结构
if collect_stable:
    stable_structs = Path(input_csv_path).parent.joinpath("stable_structs")
    if not stable_structs.exists():
        stable_structs.mkdir()

    # calypso结构预测后结构搜集
    for ent in ini_pd.stable_entries:
        print(f"look for the position of N0.{ent.entry_id} structure!")
        for ana_out_dat in Path(input_csv_path).parent.rglob("Analysis_Output.dat"):
            f = open(ana_out_dat).readlines()
            for line in f:
                patter = re.compile(r"\d+\s\(\s*%s\)" %ent.entry_id)
                result = re.search(patter, line)
                if result is not None:
                    localnumber         = re.findall("\d+", result.group())[0]
                    spacegroup_symmetry = re.search(r"\w+\(\d+\)", line).group()
                    spacegroup_number   = re.search(r"\(.*\)", spacegroup_symmetry).group().strip("()")
                    src_vaspfile = list(ana_out_dat.parent.rglob(f"UCell_{localnumber}_{spacegroup_number}.vasp"))[0]
                    dst_vaspfile = stable_structs.joinpath(f"UCell_{ent.entry_id}_{localnumber}_{spacegroup_number}.vasp")
                    shutil.copy(src_vaspfile, dst_vaspfile)
                    print(src_vaspfile)
                    break

# 搜集所有亚稳的结构
if collect_unstable:
    unstable_structs = Path(input_csv_path).parent.joinpath("unstable_structs")
    if not unstable_structs.exists():
        unstable_structs.mkdir()
    for entry in ini_entries:
        unstable_dict = {}
        energy_above_hull = ini_pd.get_e_above_hull(entry)*1000
        form_energy = ini_pd.get_form_energy(entry)
        if ebh_lower_limit <= energy_above_hull <= ebh_higher_limit:
            print(entry.entry_id)
            print(f"look for the position of N0.{entry.entry_id} structure!")
            for ana_out_dat in Path(input_csv_path).parent.rglob("Analysis_Output.dat"):
                f = open(ana_out_dat).readlines()
                for line in f:
                    patter = re.compile(r"\d+\s\(\s*%s\)" %entry.entry_id)
                    result = re.search(patter, line)
                    if result is not None:
                        localnumber         = re.findall("\d+", result.group())[0]
                        spacegroup_symmetry = re.search(r"\w+\(\d+\)", line).group()
                        spacegroup_number   = re.search(r"\(.*\)", spacegroup_symmetry).group().strip("()")
                        src_vaspfiles = list(ana_out_dat.parent.rglob(f"UCell_{localnumber}_{spacegroup_number}.vasp"))
                        if len(src_vaspfiles) == 1:
                            src_vaspfile = src_vaspfiles[0]
                        dst_vaspfile = unstable_structs.joinpath(f"UCell_{entry.entry_id}_{localnumber}_{spacegroup_number}.vasp")
                        shutil.copy(src_vaspfile, dst_vaspfile)
                        print(src_vaspfile)
                        break

usually_order="""I hope these orders will inspire you !!!
plotconvexhull_collectstructs.py -i nnconvexhull.csv -ed Ce Sr H -ebh 0.1 10 -cu -cs
"""
print(usually_order)