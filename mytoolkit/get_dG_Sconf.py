#!/usr/bin/env python

import sys
import math

print("Note: --------------------")
print("    这个脚本可以计算考虑了构型熵之基础上的能量, 即自由能, 但不是完全体的自由能, 因为这里面还没有考虑振动熵和电子熵")
print("    如果你想只取某一列或者某几列. 那么你可以在第二次运行该脚本的时候, 搭配awk来获得。比如, 你只想获得第二列(deltaH)和第三列(deltaG), 你可以这么写相关命令: ")
print("        get_dG.py | awk '{print $2, $3}'")
dH = input("请输入每原子的焓值enthalpy(单位是eV/atom)\n")
dH = float(dH)
print("dH = {:<12.8f} eV/atom".format(dH))
xA = float(input("由于这个脚本只能计算二元无序合金的构型熵, 所以你只需要输入一种元素的所占的比例xA即可, 另一种元素的比例xB=1-xA就可以获得\n"))
xB = 1 - xA
print("So percentage of B-element is {:.3f}".format(xB))
kB = 8.6173324e-5 #eV/K

Sconf_PerAtom = - kB * (xA * math.log(xA) + xB * math.log(xB))
print("So entropy is {:<12.8f} eV/K".format(Sconf_PerAtom))
print("{:<10} {:<12} {:<12}".format("T", "dH", "dG"))
for T in range(0, 3100, 100):
    dG = dH-T*Sconf_PerAtom 
    print("{:<10} {:<12.8f} {:<12.8f}".format(T, dH, dG))

