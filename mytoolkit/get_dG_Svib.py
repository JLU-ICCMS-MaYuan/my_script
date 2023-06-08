#!/usr/bin/env python

import sys
import os

from pathlib import Path
import numpy as np
import pandas as pd

def get_phonodos(all_atoms_quantity):
    """
    返回值是一个dataframe类型，是一个二维列表，第一行是freq+元素名称, 第二行
    """
    # 获得phonodos计算的输出文件
    phonon_dos_path = Path.cwd().joinpath("phonon.dos")

    # 获得元素的顺序 以及 每种元素的原子个数

    phonondos = pd.read_table(
        phonon_dos_path,
        skiprows=1,  # skiprows=1：跳过文件的第一行，即不将其作为数据的一部分进行读取。
        header=None, # header=None：不将文件的第一行作为列名，而将其视为数据。
        sep='\s+'    # sep='\s+'：使用正则表达式 \s+ 作为列之间的分隔符，表示一个或多个空格字符。
        )

    phonondos = phonondos.iloc[:, :2]
    phonondos.columns = ['freq', 'tdos']
    phonondos = phonondos[phonondos['freq'] >=0]
    
    print("检验是否将声子态密度归一化到3N")
    sumdos = phonondos['tdos'].sum()
    phonondos['tdos'] = phonondos['tdos']/sumdos*3*all_atoms_quantity
    d_freq = phonondos.loc[2, 'freq']-phonondos.loc[1, 'freq']
    print("归一化至3N前:", d_freq*sumdos)
    print("归一化至3N后:", phonondos['tdos'].sum())

    return phonondos


def get_gibbs_from_phtdos(all_atoms_quantity):

    freq_phtdos_phpdos = get_phonodos(all_atoms_quantity)

    freq_phtdos_phpdos = freq_phtdos_phpdos[freq_phtdos_phpdos["freq"]>0]  # 只保留正频率
    
    freqs = freq_phtdos_phpdos['freq'].values # 讲pandas中提取出的freqs转化为numpy类型  # freqs is circular frequency
    freqs = freqs * 2.99792458E+10  # convert unit cm-1 to Hz

    tdos   = freq_phtdos_phpdos['tdos'].values  # conserve unit is states/cm-1 = states·cm

    temperature = np.array([i for i in range(0,5100,100)])
    # temperature = np.array([0])
    # 定义积分步长 , 乘以 0.333564E-10 的目的是把积分步长变回 cm-1 
    d_freq = np.round(
        0.33356410E-10 * (freqs[-1] - freqs[0])/(len(freqs) -1), 
        decimals=6
        )  # conserve unit is cm-1
    
    # 定义计算单个声子的吉布斯自由能gibbs_i 
    def gibbs_i(freqs, T, dos):
        """\int{ G_i * g(w) } = \int{ (zpe_i + temperature_effect_i ) * g(w) dw}"""
        # h_bar = 6.582120E-16 # unit is eV·s 约化普朗克常数
        h     = 4.13566770E-15 #  unit is eV·s 普朗克常数
        k_B   = 8.61734315E-5  # unit is eV/K
        if T == 0:
            gibbs_i = (0.5 * h * freqs) * dos
        else:
            gibbs_i = (0.5 * h * freqs + k_B*T * np.log(1-np.exp(-(h*freqs)/(k_B*T)))) * dos
            # gibbs_i = k_B*T * np.log(2*np.sinh(h*freqs/(2*k_B*T))) * dos
        return gibbs_i

    T_gibbs = []
    for T in temperature:
        gibbs = gibbs_i(freqs, T, tdos) #* d_freq
        gibbs = np.sum(gibbs)
        T_gibbs.append([T, gibbs])

    with open(Path.cwd().joinpath("thermodynamics_from_phtdos.csv"), "w") as f:
        f.write("{:>5},{:>12},{:>12}\n".format("T", "gibbs(eV)", "gibbs(eV/atom)"))
        for T, gibbs in T_gibbs:
            f.write("{:>5},{:>12.8f},{:>12.8f}\n".format(T, gibbs, gibbs/all_atoms_quantity))


def get_gibbs_from_freq(all_atoms_quantity):
    """
    从声子总态密度获得声子频率然后获得自由能
        声子态密度的单位是: states/cm-1, 通过对声子态密度积分就可以得到总的振动模式数: 3N, 其中N是胞内总原子数
    """
    temperature = np.array([i for i in range(0,5100,100)])

    def gibbs_q(omegas, T, weight):
        """\int{ G_i * g(w) } = \int{ (zpe_i + temperature_effect_i ) * g(w) dw}"""
        h_bar = 6.582119514E-16 # unit is eV·s
        k_B   = 8.617343E-5  # unit is eV/K
        omegas = omegas[ omegas > 0 ]* 2*np.pi # 圆频率
        # print(len(omegas)); print(weight)
        # input(omegas)
        if T == 0:
            _gibbs_q = (1/2 * h_bar * omegas)
        else:
            _gibbs_q = (1/2 * h_bar * omegas + k_B * T * np.log(1-np.exp(-(h_bar*omegas)/(k_B*T))))
        # print(_gibbs_q); input()
        _gibbs_q = np.sum(_gibbs_q) * weight
        # print(_gibbs_q); input()
        return _gibbs_q
    dyn_paths = list(
        filter(lambda x: 'dyn' in x \
                    and 'dyn0' not in x \
                    and 'dyna2F' not in x 
                    and 'matdyn.modes' not in x \
                    and 'thermodynamics_from_dyn1.csv' not in x \
                    and "thermodynamics_from_dyns.csv" not in x \
                    and "thermodynamics_from_phtdos.csv" not in x, 
                    os.listdir(Path.cwd()))
            )
    if len(dyn_paths) == 0:
        print("There doesn't exist *.dyn*")
        sys.exit(1)

    full_dyns_paths = sorted([os.path.join(Path.cwd(), dyn_path) for dyn_path in dyn_paths])
    if not full_dyns_paths:
        print(f'WARNING: Unable to detect *dyn* files in current workpath. The program will exit')
        sys.exit(1)
    
    dyns  = [Dyn(path) for path in full_dyns_paths]

    # 从所有的q点获得频率计算吉布斯自由能
    T_gibbs_dyns = []
    for T in temperature:
        gibbs = 0
        qtot  = 0
        i = 0
        for dyn in dyns:
            qtot += dyn.weight
            gibbs_Q = gibbs_q(dyn.omegas*1.0E+12, T, dyn.weight) # dyn.omegas 的单位是THz
            gibbs += gibbs_Q
        gibbs = gibbs/qtot # 除以总的q点数
        T_gibbs_dyns.append([T, gibbs])
    with open(Path.cwd().joinpath("thermodynamics_from_dyns.csv"), "w") as f:
        f.write("{:>5},{:>12},{:>12}\n".format("T", "gibbs(eV)", "gibbs(eV/atom)"))
        for T, gibbs in T_gibbs_dyns:
            f.write("{:>5},{:>12.8f},{:>12.8f}\n".format(T, gibbs, gibbs/all_atoms_quantity))

    # 从gamma点获得频率计算吉布斯自由能
    T_gibbs_dyn1 = []
    for T in temperature:
        gibbs = 0
        gibbs_Q = gibbs_q(dyns[0].omegas*1.0E+12, T, dyns[0].weight) # dyn.omegas 的单位是THz
        gibbs += gibbs_Q
        T_gibbs_dyn1.append([T, gibbs])
    with open(Path.cwd().joinpath("thermodynamics_from_dyn1.csv"), "w") as f:
        f.write("{:>5},{:>12},{:>12}\n".format("T", "gibbs(eV)", "gibbs(eV/atom)"))
        for T, gibbs in T_gibbs_dyn1:
            f.write("{:>5},{:>12.8f},{:>12.8f}\n".format(T, gibbs, gibbs/all_atoms_quantity))


class Dyn(object):
    """
    Parses a single *dyn*.elph* file,
    returning all it contains:
    q-point, lambdas, gammas and squared frequencies
    """
    lines = list()
    q_point = tuple()
    weight = float()


    def __init__(self, path):
        with open(path) as read_obj:
            lines = read_obj.readlines()
            read_obj.close()
        self.lines = lines
        # print(f'q = ({", ".join(["%.3f" % round(_q, 3) for _q in self.q_point])}) '
            #   f'with number of q in the star {int(self.weight)}')
        
    @property
    def q_point(self):
        q_idx = int()
        for idx, line in enumerate(self.lines):
            if 'q = (' in line:
                q_idx = idx
                break
        self._q_point = tuple(float(x) for x in self.lines[q_idx].split()[-4:-1])
        return self._q_point
    
    @property
    def weight(self):
        self._weight = 0
        for line in self.lines:
            if 'Diagonalizing the dynamical matrix' in line: # 这一行前面都是等价的q点坐标，这一行后面是针对这些等价的坐标选取一个进行动力学矩阵对角化
                break
            if 'q = (' in line:
                self._weight = self._weight + 1
        # self.weight = float(self.lines[2].split()[1])
        return self._weight

    @property
    def omegas(self):
        self._omegas = []
        for line in self.lines:
            if 'freq (' in line:
                omega = float(line.split()[4])
                self._omegas.append(omega)
        self._omegas = np.array(self._omegas)
        return self._omegas


if __name__ == "__main__":
    print("make sure your file name is phonon.dos")
    atoms_amount = float(input("you need to input number of atoms in a cell\n"))

    try:
        get_gibbs_from_phtdos(atoms_amount)
        print("the results was save in thermodynamics_from_phtdos.csv")
    except:
        print("get free-energy from phonon.dos failed!")

    try:
        get_gibbs_from_freq(atoms_amount)
        print("the results was save in thermodynamics_from_dyns.csv and thermodynamics_from_dyn1.csv")
    except:
        print("get free-energy from *.dyn* failed! Make sure in current path, *dyn* exists! ")
