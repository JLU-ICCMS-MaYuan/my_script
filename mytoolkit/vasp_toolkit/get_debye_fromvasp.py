#!/usr/bin/env python

import sys
import os

from pathlib import Path
import numpy as np
import pandas as pd


def get_phonodos(all_atoms_quantity):
    """
    返回值是一个dataframe类型, 是一个二维列表, 第一行是freq+元素名称, 第二行
    """
    # 获得phonodos计算的输出文件
    phonon_dos_path = Path.cwd().joinpath("phonon.dos")

    # 获得元素的顺序 以及 每种元素的原子个数

    phonondos = pd.read_table(
        phonon_dos_path,
        skiprows=1,  # skiprows=1：跳过文件的第一行, 即不将其作为数据的一部分进行读取。
        header=None, # header=None：不将文件的第一行作为列名, 而将其视为数据。
        sep='\s+'    # sep='\s+'：使用正则表达式 \s+ 作为列之间的分隔符, 表示一个或多个空格字符。
        )

    phonondos = phonondos.iloc[:, :2]
    phonondos.columns = ['freq', 'tdos']
    phonondos = phonondos[phonondos['freq'] >=0]

    # print("检验是否将声子态密度归一化到3N")
    # sumdos = phonondos['tdos'].sum()
    # phonondos['tdos'] = phonondos['tdos']/sumdos*3*all_atoms_quantity
    # d_freq = phonondos.loc[2, 'freq']-phonondos.loc[1, 'freq']
    # print("归一化至3N前:", d_freq*sumdos)
    # print("归一化至3N后:", phonondos['tdos'].sum())

    return phonondos


def set_Debye_frequency(phonondos, num_atoms, freq_max_fit=None):
    """Calculate a kind of Debye frequency."""
    try:
        from scipy.optimize import curve_fit
    except ImportError:
        print("You need to install python-scipy.")
        sys.exit(1)

    def Debye_dos(freq, a):
        return a * freq**2

    freq = phonondos['freq']
    tdos = phonondos['tdos']

    freq_min = freq.min()
    freq_max = freq.max()
    print("freq_min", freq_min)
    print("freq_max", freq_max)
    if freq_max_fit is None:
        N_fit = int(len(freq) / 4.0)  # Hard coded
    else:
        N_fit = int(
            freq_max_fit / (freq_max - freq_min) * len(freq)
        )
    popt, pcov = curve_fit(Debye_dos, freq[0:N_fit], tdos[0:N_fit])
    # popt 是最佳拟合参数的估计值
    # pcov 是参数估计的协方差矩阵
    a2 = popt[0]
    _freq_Debye = (3 * 3 * num_atoms / a2) ** (1.0 / 3)
    _Debye_fit_coef = a2

    cm_1toThz = 0.0299792458
    cm_1toKelvin = 1.438776877
    ThzToKelvin = 47.99243073
    print("_Debye_fit_coef", _Debye_fit_coef)
    print("_freq_Debye",_freq_Debye)
    print("_temperture_Debye", _freq_Debye*ThzToKelvin)

if __name__ == "__main__":


    print("Note: --------------------")
    print("   get_debye_fromvasp.py [dosfile] [atom_num_in_a_cell] [max_frequency]")
    print("   The script does not read the first line of the [dosfile], \n\
              only reads from the second line to the last line.\n \
              The first column and the second column are read, \n\
              where the first column must be the frequency in units of THz, \n\
              and the second column must be the total phonon density of states.\n")
    dosfile = sys.argv[1]
    num_atoms = float(sys.argv[2])
    try:
        freq_max_fit = float(sys.argv[3])
    except:
        freq_max_fit = None
    phonondos = get_phonodos(num_atoms)
    set_Debye_frequency(phonondos, num_atoms, freq_max_fit)
    