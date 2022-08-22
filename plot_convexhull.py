#!/public/home/mayuan/miniconda3/envs/cage/bin/python3
#!/work/home/mayuan/miniconda3/envs/cage/bin/python3

'''
展示画图结果, 可动态调节:
    plot_convexhull.py 指定csv文件路径 show
    plot_convexhull.py ./B-H.csv show
保存画图结果,
    plot_convexhull.py 指定csv文件路径 save 
    plot_convexhull.py ./B-H.csv save
将落在convex hull上的稳定结构拿出来放在目录stable_structs中
    plot_convexhull.py 指定csv文件路径 collect
    plot_convexhull.py nnconvexhull.csv collect
'''

import sys
import os
import re
import shutil
from pathlib import Path

import pandas as pd
from pymatgen.analysis.phase_diagram import PDEntry, PhaseDiagram, PDPlotter
from pymatgen.core.composition import Composition
from pymatgen.core.periodic_table import DummySpecies, Element

input_csv = sys.argv[1]; input_csv_path = Path(input_csv)
convexhull_data = pd.read_csv(input_csv_path, header=0, sep=',')
if "save" in sys.argv:
    save_pnd = True
else:
    save_pnd = False
if "show" in sys.argv:
    show_pnd = True
else:
    show_pnd = False
if "collect" in sys.argv:
    collect_stable = True
else:
    collect_stable = False

ini_entries = []
for idx, row in convexhull_data.iterrows():
    comp = Composition(row['formula'])
    enth = row['enthalpy']
    entry_id = row['Number']
    _entry = PDEntry(comp, enth)
    _entry.entry_id = entry_id
    ini_entries.append(_entry)

ini_pd = PhaseDiagram(ini_entries)
print(f" reference material {ini_pd.el_refs}\n");
#获得落在convex hull上的结构的化学式配比, 编号(编号用来索引搜索到的结构)
stable_structs_amount = 0
for ent in ini_pd.stable_entries:
    form_energy = ini_pd.get_form_energy(ent)
    print("stoichiometry: {:<10}  No.{:<10}, its form energy is {:<10}".format(
        ent.name, 
        ent.entry_id,
        form_energy,
        ))
    stable_structs_amount += 1


print(f"stable structures on the convex hull is {stable_structs_amount - len(ini_pd.el_refs)}\n")

for ent in ini_entries:
    ehull       = ini_pd.get_e_above_hull(ent)*1000
    form_energy = ini_pd.get_form_energy(ent)
    if 0.0< ehull <= 1:
        print("stoichiometry: {:<10}  No.{:<10} is above convell hull {:<10}".format(
            ent.name, 
            ent.entry_id,
            ehull,
            ))

if save_pnd:
    plotter = PDPlotter(ini_pd, show_unstable=0.2, backend='matplotlib')
    plotter.write_image('pd.png', image_format='png')

if show_pnd:
    plotter = PDPlotter(ini_pd, show_unstable=0.2, backend='plotly')
    plotter.show()

print("\n")

if collect_stable:
    stable_structs = Path(input_csv_path.parent).joinpath("stable_structs")
    if not stable_structs.exists():
        stable_structs.mkdir()

    for ent in ini_pd.stable_entries:
        print(f"look for the position of N0.{ent.entry_id} structure!")
        for ana_out_dat in input_csv_path.parent.rglob("Analysis_Output.dat"):
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





