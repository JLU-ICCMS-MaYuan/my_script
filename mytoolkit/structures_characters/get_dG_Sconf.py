#!/usr/bin/env python

import sys
import math

from fractions import Fraction

def multinary_alloy(dG_file, natoms, ratios):

    print("Ratio in every elements:")
    for rt in ratios:
        print("{:<6.5f}".format(rt))
    if sum(ratios) - 1.0 < 0.001:
        print("Good input ratio, summary result is 1")
    else:
        print("比例不对, 和{} != 1".format(sum(ratios)))
        sys.exit()

    kB = 8.6173324e-5 #eV/K

    Sconf_Performula = 0
    for rt in ratios:
        Sconf_Performula += -rt * math.log(rt)
    Sconf_Peratom   = Sconf_Performula/natoms

    print("Sconf = {:<12.8f} kB/formula = {:<12.8f} kB/atom".format(Sconf_Performula, Sconf_Peratom))
    print("Sconf = {:<12.8f} eV/K/formula = {:<12.8f} eV/K/atom".format(Sconf_Performula*kB, Sconf_Peratom*kB))
    print("Sconf = {:<12.8f} meV/K/formula = {:<12.8f} meV/K/atom".format(Sconf_Performula*kB*1000, Sconf_Peratom*kB*1000))

    print("{:<10} {:<12}".format("T", "-TS (eV/atom)"), file=dG_file)
    for T in range(0, 3100, 100):
        minus_ts_item = -T*Sconf_Peratom*kB
        print("{:<10} {:>12.8f}".format(T, minus_ts_item), file=dG_file)
    

if __name__ == "__main__":
    print("Note: --------------------")
    print("    这个脚本可以计算考虑了构型熵之基础上的能量, 即自由能, 但不是完全体的自由能, 因为这里面还没有考虑振动熵和电子熵")
    print("    如果你想只取某一列或者某几列. 那么你可以在第二次运行该脚本的时候, 搭配awk来获得")
    print("        get_dG.py | awk '{print $1}'")

    info1 = """请输入总原子数
例如: LaBeH8 在La的占位上合金化, LaBeH8原胞的总原子数为10, natom=10;
      hcp-Fe 在Fe的占位上合金化, hcp-Fe原胞的总原子数为2, 但是其分子式中原胞数为1, 所以natom=1;
"""

    natoms = int(input("请输入总原子数\n"))

    info2 = """请输入你要对几个占位进行合金化,每个占位混入多少种金属.
例如: LaBeH8 仅在La的占位上合金化, 用4种元素混合, 则输入每种元素的混合百分比: 0.25 0.25 0.25 0.25
      hcp-Fe  在两个Fe的占位上合金化, 用5种元素混合, 则输入每种元素的混合百分比: 6/10 1/10 1/10 1/10, 在指定百分比时并没有指明哪个元素占据哪个占位,全部都是随机的占据
      保证你输入的百分比加和为1, 否在程序报错
      """

    ratios = input(info2).split()
    float_ratios = []
    for rt in ratios:
        try:
            value = float(Fraction(rt))  # 将分数字符串转换为 Fraction 对象，然后转换为浮点数
        except ValueError as e:
            value = float(rt)
        float_ratios.append(value)
    dG_file = open("dG.dat", 'w')
    print("natom={:<12}".format(natoms), file=dG_file)
    multinary_alloy(dG_file, natoms, float_ratios)
