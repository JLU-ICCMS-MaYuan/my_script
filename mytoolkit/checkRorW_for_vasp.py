#!/usr/bin/env python3
import os
import re

print("当前工作路径:")
print(os.getcwd())


print("----检查赝势----")
title = os.popen("grep TITEL POTCAR").read().split('\n')
title = [t for t in title if t]
zval = os.popen("grep ZVAL POTCAR").read().split('\n')
zval = [z for z in zval if z]
for t, z in zip(title, zval):
    zv = z.split()[5]
    print("{:<50}    ZVAL={:<10}".format(t, zv))

print("----POSCAR元素顺序----")
element = os.popen("sed -n '6p' POSCAR").read().split()
for e in element:
    print("{:>5}".format(e))

print("----检查INCAR----")
press = os.popen("grep PSTRESS INCAR").read()
print("press    = {} GPa".format(int(re.search(r'\d+', press).group())/10))
os.system("grep ISMEAR INCAR")
os.system("grep SIGMA INCAR")
os.system("grep ENCUT INCAR")
os.system("grep KSPACING INCAR")

print("----检查POSCAR和CONTCAR的对称性----")
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
struct = Structure.from_file("POSCAR")
spgaly = SpacegroupAnalyzer(struct, symprec = 1E-5)
symbol = spgaly.get_space_group_symbol()
spgnum = spgaly.get_space_group_number()
print("POSCAR: {} ({})".format(symbol, spgnum))

struct = Structure.from_file("CONTCAR")
spgaly = SpacegroupAnalyzer(struct, symprec = 1E-5)
symbol = spgaly.get_space_group_symbol()
spgnum = spgaly.get_space_group_number()
print("CONTCAR: {} ({})".format(symbol, spgnum))