#!/usr/bin/env python

import sys
import math

print("Note: --------------------")
print("    这个脚本可以计算考虑了构型熵之基础上的能量, 即自由能, 但不是完全体的自由能, 因为这里面还没有考虑振动熵和电子熵")
xA = float(input("由于这个脚本只能计算二元无序合金的构型熵, 所以你只需要输入一种元素的所占的比例xA即可, 另一种元素的比例xB=1-xA就可以获得\n"))
xB = 1 - xA
print("So percentage of B-element is {:.3f}".format(xB))
kB = 8.6173324e-5 #eV/K

Sconf_PerAtom = - kB * (xA * math.log(xA) + xB * math.log(xB))
print("So entropy is {:<12.8f} eV/K".format(Sconf_PerAtom))
print("{:<10} {:<12}".format("T", "T*S"))
for T in range(0, 3100, 100):
    ts_iterm = T*Sconf_PerAtom 
    print("{:<10} {:<12.8f}".format(T, ts_iterm))

dH_formation = float(input("输入delta_H(反应物)-delta(生成物)/N的值, 即输入形成焓(千万注意你输入的必须是:eV/atom !!!!!!!)\n"))
if dH_formation:
    print("{:<10} {:<12} {:<12}".format("T", "dH", "dG"))
    for T in range(0, 3100, 100):
        ts_iterm = T*Sconf_PerAtom 
        dG = dH_formation+ts_iterm
        print("{:<10} {:<12.8f} {:<12.8f}".format(T, dH_formation, dG))

