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
from pprint import pprint
from argparse import ArgumentParser, RawTextHelpFormatter

import pandas as pd
from pymatgen.analysis.phase_diagram import PDEntry, PDPlotter, CompoundPhaseDiagram
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
    "-ed",
    "--endnotes",
    action="store",
    default=None,
    dest="endnotes",
    nargs="+",
    help="如果使用该参数, 需要用户指定某几个化合物作为端点值"
)
parser.add_argument(
    "-dei",
    "--destination-entry-index",
    action="store",
    default=[-1],
    nargs="+",
    dest="dst_entry_index",
    type=int,
    help="该参数用于输出指定entry的分解路径, 使用者需要事先了解相应化合物的entry编号\n"
        "一般情况看nnconvexhull.csv文件中目标化合物在第几行就可以, 例如:\n"
        "对于最后一行的化合物可以用-1来代表最后一个entry\n"
)
parser.add_argument(
    "-save",
    "--save-png",
    action="store",
    default=0.0,
    type=float,
    dest="save_png",
    help="是否保存图片, 数值为显示ebh的结构的上限, 这里保存的图片是pymatgen自动生成的。\n"
         "保存convexhull图片到当前执行命令的路径下\n"
)
parser.add_argument(
    "-show",
    "--show-png",
    action="store_true",
    default=0.0,
    type=float,
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
    type=float,
    dest="EnthalpyAboveHullValue",
    default=[0, 10.0],
    nargs="+",
    help="高于convex hull xxx emV 的能量的上限\n"
        "在 ebh_lower_limit ~ ebh_higher_limit 这个范围内的亚稳的结构确定出来\n"
        "所以说要输入两个值，第一个值是下陷，第二个值是上限。\n"
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
ebh_lower_limit = args.EnthalpyAboveHullValue[0]
ebh_higher_limit = args.EnthalpyAboveHullValue[1]
endnotes = args.endnotes
dst_entry_index = args.dst_entry_index
# 生成 凸包图对象
# convexhull_data = pd.read_csv(input_csv_path, header=0, sep=',') #  header表示第一行为标题行
print("读入文件中的能量和化学式 (注意：能量必须是化学式的能量，不是每原子的能量) ")
try:
    # 该情况处理的文件：
    # Number  formula     enthalpy
    # 1       Ax1By1Cz1   -1.34343
    # 2       Ax2By2C2z   -2.324324
    convexhull_data = pd.read_table(input_csv_path, header=0, sep=',') #  header表示第一行为标题行 sep='\s+'表示使用多个空格作为分隔符
    ini_entries = []
    for idx, row in convexhull_data.iterrows():
        comp = Composition(row['formula'])
        num_at = comp.num_atoms
        enth = row['enthalpy']*num_at
        entry_id = row['Number']
        _entry = PDEntry(comp, enth, attribute={"entry_id":entry_id})
        #print(_entry.attribute['entry_id'])
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
        _entry = PDEntry(comp, enth, attribute={"entry_id":entry_id})
        #print(_entry.attribute['entry_id'])
        ini_entries.append(_entry)


# 建立相图
terminal_comps = [Composition(formula) for formula in endnotes]
ini_pd = CompoundPhaseDiagram(
    entries = ini_entries,
    terminal_compositions = terminal_comps,
    normalize_terminal_compositions = True,
    )
    
# 输出参考单质或化合物
print(ini_pd.terminal_compositions)
for re in ini_pd.el_refs:
    print(re.symbol)
print(f"reference endnotes\n")


# 获得新的变换坐标后的entry
# trans_entries = []
# for entry in ini_entries:
#     new_entries, sp_mapping = TransformedPDEntry(entry, sp_mapping)(entries=entry, terminal_compositions=terminal_comps)
#     trans_entries.append(new_entries)

# 获得 落在convex hull上的稳定结构的 化学式配比, 编号, 形成焓(编号用来索引搜索到的结构)
# 获得 落在convex hull 上稳定结构的 csv 文件
stable_list = []
stable_structs_amount = 0
for entry in ini_pd.stable_entries:
    stable_dict = {}
    form_energy = ini_pd.get_form_energy(entry)
    stable_dict["Number"] = entry.attribute['entry_id']
    stable_dict["dummyformula"] = entry.composition.formula.replace(' ', '')
    stable_dict["formula"] = entry.name
    stable_dict["enthalpy"] = 0.0
    stable_list.append(stable_dict)
    stable_structs_amount += 1
print(f"stable structures on the convex hull is {stable_structs_amount - len(ini_pd.el_refs)}\n")
print("{:<10}  {:<10}  {:<20}   {:<10}".format("Number", "formula", "dummyformula", "enthalpy(eV/atom)"))
sorted_dict_list = sorted(stable_list, key=lambda d: d['formula'])
for unstable_dict in sorted_dict_list:
    print("{:<10}  {:<10}  {:<20}  {:<10.4f}".format(
        unstable_dict["Number"],
        unstable_dict["formula"], 
        unstable_dict["dummyformula"], 
        unstable_dict["enthalpy"],
        ))
stable_pd = pd.DataFrame(stable_list)
stable_pd.to_csv("stable.csv", index=False)


# 获得 高于convex hull 上亚稳结构的 csv 文件
unstable_list = []
unstable_structs_amount = 0
for entry in ini_pd.unstable_entries:
    unstable_dict = {}
    energy_above_hull = ini_pd.get_e_above_hull(entry)*1000
    form_energy = ini_pd.get_form_energy(entry)
    if ebh_lower_limit < energy_above_hull <= ebh_higher_limit: # 这里取高于convex hull 能量在0~50个meV范围内的亚稳结构
        unstable_dict["Number"]  = entry.attribute['entry_id']
        unstable_dict["dummyformula"] = entry.composition.formula.replace(' ', '')
        unstable_dict["formula"] = entry.name
        #!!!!!!!!!!!!!!!!!特别注意这里的单位是meV!!!!!!!!!!!!!!!!!!!!!!!
        unstable_dict["enthalpy"] = ini_pd.get_e_above_hull(entry) 
        unstable_list.append(unstable_dict)
        unstable_structs_amount += 1
sorted_dict_list = sorted(unstable_list, key=lambda d: d['enthalpy'])
print(f"unstable structures above the convex hull {ebh_lower_limit}-{ebh_higher_limit} meV is {unstable_structs_amount}\n")
print("{:<10}  {:<10}  {:<20}  {:<10}".format("Number", "formula", "dummyformula", "enthalpy(eV/atom)"))
for unstable_dict in sorted_dict_list:
    print("{:<10}  {:<10}  {:<20}  {:<10.4f}".format(
        unstable_dict["Number"],
        unstable_dict["formula"], 
        unstable_dict["dummyformula"], 
        unstable_dict["enthalpy"],
        ))
unstable_pd = pd.DataFrame(unstable_list)
unstable_pd.to_csv("unstable.csv", index=False)



plotter = PDPlotter(ini_pd, show_unstable=save_pnd, backend='matplotlib')
plotter.write_image('pd_cpd.png', image_format='png')


if dst_entry_index:
    for dei in dst_entry_index:
        dei = int(dei)
        dst_entry = ini_pd.entries[dei]
        print(f"\n------------dst_entry[{-1}]={dst_entry.original_entry.composition.reduced_formula}------------")
    
        decomp_path = ini_pd.get_decomp_and_e_above_hull(dst_entry)[0]
        decomp_eabovehull = ini_pd.get_decomp_and_e_above_hull(dst_entry)[1]
        form_energy = ini_pd.get_form_energy_per_atom(dst_entry)
        print("{:<15} {:<15} {:<20} {:<20} {:<20}".format('formula_trans', 'formula', 'energy_per_atom', 'energy', 'ratio'))
        for key, value in decomp_path.items():
            print("{:<15} {:<15} {:<20} {:<20} {:<20}".format(key.composition.reduced_formula, key.original_entry.composition.reduced_formula, key.energy_per_atom, key.energy, value))
        print("   e_above_hull = {:>10.6f} ev/atom".format(decomp_eabovehull))
        print("formed_enthalpy = {:>10.6f} ev/atom".format(form_energy))
