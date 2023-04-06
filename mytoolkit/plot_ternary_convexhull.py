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
    "-save",
    "--save-png",
    action="store_true",
    default=False,
    dest="save_png",
    help="是否保存图片, 这里保存的图片是pymatgen自动生成的。\n"
        "保存convexhull图片到当前执行命令的路径下\n"
)
parser.add_argument(
    "-show",
    "--show-png",
    action="store_true",
    default=False,
    dest="show_png",
    help="是否展示图片, 以网页的形式展示出一个可以动态调节详细程度的可交互图片\n",
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
    type=int,
    dest="EnthalpyAboveHullValue",
    help="高于convex hull xxx emV 的能量的上限\n"
        "在0 ~ EnthalpyAboveHullValue 这个范围内的亚稳的结构确定出来\n"
)
parser.add_argument(
    "-hand",
    "-hand-plot-dat",
    action="store",
    default=None,
    dest="hand_plot_dat",
    nargs="+",
    help="导出精细的数据方便origin绘制"
        "如果是一个三元的体系, 该参数应该这样设置：\n"
        "   -hand Xn Ym Zq  \n"
        "########### 这里一定要注意 #############\n"
        "如果单质X 有 n 个原子,一定要写对 n 的个数\n"
        "#######################################\n"
        "然后将会输出3个文件: stable.csv, unstable.csv, for_origin_plot.csv\n"
        "stable.csv 存储着所有稳定结构的焓值\n"
        "unstable.csv 存储着所有压稳结构高于convex hull的能量值\n"
        "for_origin_plot.csv 存储着所有的原始数据"
)
args = parser.parse_args()

input_csv_path = args.input_csv_path
save_pnd = args.save_png
show_pnd = args.show_png
collect_stable = args.collect_stable
collect_unstable = args.collect_unstable
hand_plot_dat = args.hand_plot_dat
EnthalpyAboveHullValue = args.EnthalpyAboveHullValue
# 生成 凸包图对象
# convexhull_data = pd.read_csv(input_csv_path, header=0, sep=',') #  header表示第一行为标题行
convexhull_data = pd.read_table(input_csv_path, header=0, sep='\s+') #  header表示第一行为标题行
ini_entries = []
for idx, row in convexhull_data.iterrows():
    comp = Composition(row['formula'])
    enth = row['enthalpy']
    entry_id = row['Number']
    _entry = PDEntry(comp, enth)
    _entry.entry_id = entry_id
    ini_entries.append(_entry)
ini_pd = PhaseDiagram(ini_entries)

# 输出参考单质或化合物
print(f" reference material {ini_pd.el_refs}\n");



# 获得 落在convex hull上的稳定结构的 化学式配比, 编号, 形成焓(编号用来索引搜索到的结构)
# 获得 落在convex hull 上稳定结构的 csv 文件
stable_list = []
stable_structs_amount = 0
for entry in ini_pd.stable_entries:
    stable_dict = {}
    form_energy = ini_pd.get_form_energy(entry)
    stable_dict["Number"] = entry.entry_id
    stable_dict["formula"] = entry.composition.formula
    print(entry.composition.formula)
    stable_dict["enthalpy"] = 0.0
    stable_list.append(stable_dict)
    stable_structs_amount += 1
print(f"stable structures on the convex hull is {stable_structs_amount - len(ini_pd.el_refs)}\n")
stable_pd = pd.DataFrame(stable_list)
stable_pd.to_csv("stable.csv", index=False)

# 获得 高于convex hull  0~50mev 上亚稳结构的 csv 文件
unstable_list = []
unstable_structs_amount = 0
for entry in ini_entries:
    unstable_dict = {}
    energy_above_hull = ini_pd.get_e_above_hull(entry)*1000
    form_energy = ini_pd.get_form_energy(entry)
    if 0.0 < energy_above_hull <= EnthalpyAboveHullValue: # 这里取高于convex hull 能量在0~50个meV范围内的亚稳结构
        print("stoichiometry: {:<10}  No.{:<10} is above convell hull {:<10}".format(
            entry.name, 
            entry.entry_id,
            energy_above_hull,"meV",
            ))
        unstable_dict["Number"]  = entry.entry_id
        unstable_dict["formula"] = entry.composition.formula
        #!!!!!!!!!!!!!!!!!特别注意这里的单位是meV!!!!!!!!!!!!!!!!!!!!!!!
        unstable_dict["enthalpy"] = ini_pd.get_e_above_hull(entry) 
        unstable_list.append(unstable_dict)
        unstable_structs_amount += 1
print(f"unstable structures above the convex hull 0-50 meV is {unstable_structs_amount}\n")
unstable_pd = pd.DataFrame(unstable_list)
unstable_pd.to_csv("unstable.csv", index=False)



if save_pnd:
    plotter = PDPlotter(ini_pd, show_unstable=EnthalpyAboveHullValue*0.001, backend='matplotlib')
    plotter.write_image('pd.png', image_format='png')

if show_pnd:
    plotter = PDPlotter(ini_pd, show_unstable=EnthalpyAboveHullValue*0.001, backend='plotly')
    plotter.show()

print("\n")

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
    # 自定义结构优化后结构搜集
    for ent in ini_pd.stable_entries:
        src_vaspfile = Path.cwd().joinpath(str(ent.entry_id), "POSCAR")
        dst_vaspfile = stable_structs.joinpath(f"{ent.entry_id}_{ent.composition.formula.replace(' ', '')}_.vasp")
        shutil.copy(src_vaspfile, dst_vaspfile)
        print(src_vaspfile)

# 搜集所有亚稳的结构
if collect_unstable:
    unstable_structs = Path(input_csv_path).parent.joinpath("unstable_structs")
    if not unstable_structs.exists():
        unstable_structs.mkdir()
    for entry in ini_entries:
        unstable_dict = {}
        energy_above_hull = ini_pd.get_e_above_hull(entry)*1000
        form_energy = ini_pd.get_form_energy(entry)
        if 0.0 < energy_above_hull <= 50.0:
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
                        src_vaspfile = list(ana_out_dat.parent.rglob(f"UCell_{localnumber}_{spacegroup_number}.vasp"))[0]
                        dst_vaspfile = unstable_structs.joinpath(f"UCell_{entry.entry_id}_{localnumber}_{spacegroup_number}.vasp")
                        shutil.copy(src_vaspfile, dst_vaspfile)
                        print(src_vaspfile)
                        break

    # 自定义结构优化后结构搜集
    for ent in ini_pd.ini_entries:
        unstable_dict = {}
        energy_above_hull = ini_pd.get_e_above_hull(entry)*1000
        form_energy = ini_pd.get_form_energy(entry)        
        if 0.0 < energy_above_hull <= 50.0:
            src_vaspfile = Path.cwd().joinpath(str(ent.entry_id), "POSCAR")
            dst_vaspfile = unstable_structs.joinpath(f"{ent.entry_id}_{ent.composition.formula.replace(' ', '')}_.vasp")
            shutil.copy(src_vaspfile, dst_vaspfile)
            print(src_vaspfile)


if hand_plot_dat:

    for ele in hand_plot_dat:
        stable_pd.insert(2, ele, 0)
    for idx, row in stable_pd.iterrows():
        entry_id = row['Number']
        comp = Composition(row['formula']) 
        for ele in comp.elements:
            ele_name = ele.name
            atoms_num= comp.to_data_dict['unit_cell_composition'][ele_name]
            stable_pd.loc[idx, ele_name] = atoms_num
    stable_pd.to_csv("stable.csv", index=False)


    for ele in hand_plot_dat:
        unstable_pd.insert(2, ele, 0)
    for idx, row in unstable_pd.iterrows():
        entry_id = row['Number']
        comp = Composition(row['formula']) 
        for ele in comp.elements:
            ele_name = ele.name
            atoms_num= comp.to_data_dict['unit_cell_composition'][ele_name]
            unstable_pd.loc[idx, ele_name] = atoms_num
    unstable_pd.to_csv("unstable.csv", index=False)

    total_pd = pd.concat((stable_pd, unstable_pd), axis=0)
    #total_pd = total_pd.drop(columns=['Number', 'formula']) # 删除列
    total_pd.to_csv("for_origin_plot.csv", index=False)





