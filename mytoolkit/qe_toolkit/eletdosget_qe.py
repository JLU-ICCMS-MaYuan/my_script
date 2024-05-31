#!/usr/bin/env python3
import os
import sys
import numpy as np
import pandas as pd

file = sys.argv[1]

print("NOTE --------------------")
print("    This scrpit can get N(Ef) from *.tdos file obtained by nscf-calculation")
print("    You need to specify the *.tdos filename, it can be known from the parameter of `fildos = 'La1Ce1Be2H16.tdos',` in `eletdos.in` file \n")


efermi = os.popen("grep Fermi scffit.out | tail -n 1 | awk '{print $5}'").read()
print("NOTE --------------------")
print("    Get fermi-energy from scffit.out")
print(f"    efermi={efermi}")

# 载入tdos
data = np.loadtxt(file, comments='#')
# 计算data第一列中每个数与费米能级efermi的差值的绝对值
data[:,0] = data[:,0] - float(efermi)

idx = np.abs(data[:, 0]).argmin()
dst_data = data[idx, :]

print("NOTE --------------------")
print("    N(Ef) = {:.6f} states/eV/(Unit Cell) from nscf-calculation\n          = {:.6f} states/spin/Ry/(Unit Cell) after unit convertion".format(dst_data[1], dst_data[1]/2/0.0734986443513))
print("    The unit of N(Ef) in lambda.out is `states/spin/Ry/(Unit Cell)`")

print("NOTE --------------------")
print("    The eletdos_fromqe.csv can be used in origin to plot-picture, and the efermi has been reduced !!!!!!")
pd.DataFrame(data, columns=['energy(eV)', 'dos(states/eV/unitcell)', 'dos_int']).to_csv(
    "eletdos_fromqe.csv",
    index=False, #  index=False表示不写入行索引
)