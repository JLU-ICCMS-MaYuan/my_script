#!/usr/bin/env python

import sys
import math



def binary_alloy(dH, dG_file):
    print("计算2元无序合金的构型熵, 所以你只需要输入一种元素的所占的比例xA即可, 另一种元素的比例xB=1-xA就可以获得\n")
    print("输入的数精确到小数点后四位")
    xA = float(input("xA="))
    xB = 1 - xA
    print("xA={:<12.8f}, xB={:<12.8f}".format(xA, xB))
    kB = 8.6173324e-5 #eV/K

    Sconf_PerAtom = - kB * (xA * math.log(xA) + xB * math.log(xB))
    print("So entropy is {:<12.8f} eV/K".format(Sconf_PerAtom))
    print("{:<10} {:<12} {:<12} {:<12}".format("T", "dH", "TS", "dG"), file=dG_file)
    for T in range(0, 3100, 50):
        ts_item = T*Sconf_PerAtom 
        dG = dH-ts_item
        print("{:<10} {:<12.8f} {:<12.8f} {:<12.8f}".format(T, dH, ts_item, dG), file=dG_file)

def ternary_alloy(dH, dG_file):
    print("计算3元无序合金的构型熵, 所以你只需要输入2种元素的所占的比例xA, xB即可, 另一种元素的比例xC=1-xA-xB就可以获得")
    print("输入的数精确到小数点后四位")
    xA = float(input("xA="))
    xB = float(input("xA="))
    xC = 1 - xA - xB
    print("xA={:<12.8f}, xB={:<12.8f}, xC={:<12.8f}".format(xA, xB, xC), file=dG_file)
    if xA+xB+xC - 1.0 < 0.001:
        print("xA+xB+xC+xD=1")
    else:
        print("xA+xB+xC+xD={}".format(xA+xB+xC), file=dG_file)
        sys.exit(1)

    kB = 8.6173324e-5 #eV/K
    Sconf_PerAtom = - kB * (xA * math.log(xA) + xB * math.log(xB) + xC * math.log(xC))
    print("So entropy is {:<12.8f} eV/K".format(Sconf_PerAtom))
    print("{:<10} {:<12} {:<12} {:<12}".format("T", "dH", "TS", "dG"), file=dG_file)
    for T in range(0, 3100, 50):
        ts_item = T*Sconf_PerAtom 
        dG = dH-ts_item
        print("{:<10} {:<12.8f} {:<12.8f} {:<12.8f}".format(T, dH, ts_item, dG), file=dG_file)

def quaternary_alloy(dH, dG_file):
    print("计算4元无序合金的构型熵, 所以你只需要输入一种元素的所占的比例xA, xB, xC即可, 另一种元素的比例xD=1-xA-xB-xC就可以获得\n")
    print("输入的数精确到小数点后四位")
    xA = float(input("xA="))
    xB = float(input("xB="))
    xC = float(input("xC="))
    xD = 1 - xA - xB -xC

    print("xA={:<12.8f}, xB={:<12.8f}, xC={:<12.8f}, xD={:<12.8f}".format(xA, xB, xC, xD))
    if xA+xB+xC+xD - 1.0 < 0.001:
        print("xA+xB+xC+xD=1")
    else:
        print("xA+xB+xC+xD={}".format(xA+xB+xC+xD))
        sys.exit(1)

    kB = 8.6173324e-5 #eV/K

    Sconf_PerAtom = - kB * (xA * math.log(xA) + xB * math.log(xB) + xC * math.log(xC) + xD * math.log(xD))
    print("So entropy is {:<12.8f} eV/K".format(Sconf_PerAtom))
    print("{:<10} {:<12} {:<12} {:<12}".format("T", "dH", "TS", "dG"), file=dG_file)
    for T in range(0, 3100, 50):
        ts_item = T*Sconf_PerAtom 
        dG = dH-ts_item
        print("{:<10} {:<12.8f} {:<12.8f} {:<12.8f}".format(T, dH, ts_item, dG), file=dG_file)
    

if __name__ == "__main__":
    print("Note: --------------------")
    print("    这个脚本可以计算考虑了构型熵之基础上的能量, 即自由能, 但不是完全体的自由能, 因为这里面还没有考虑振动熵和电子熵")
    print("    如果你想只取某一列或者某几列. 那么你可以在第二次运行该脚本的时候, 搭配awk来获得。比如, 你只想获得第二列(deltaH)和第三列(deltaG), 你可以这么写相关命令: ")
    print("        get_dG.py | awk '{print $2, $3}'")
    dH = input("请输入每原子的焓值enthalpy(单位是eV/atom)\n")
    dH = float(dH)
    print("dH = {:<12.8f} eV/atom".format(dH))

    yuan = int(input("请输入你要计算的是几元合金(输入2,3,4数字即可):\n"))
    dG_file = open("dG.dat", 'w')
    if yuan == 2:
        binary_alloy(dH, dG_file)
    elif yuan == 3:
        ternary_alloy(dH, dG_file)
    elif yuan == 4:
        quaternary_alloy(dH, dG_file)