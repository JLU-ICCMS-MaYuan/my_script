#!/usr/bin/env python
# _*- coding: utf-8 -*-
# Author: Li Zhu < zhulipresent@gmail.com>
# modified by Hanyu Liu <lhy@calypso.cn> at Dec 13th, 2014
# modified by Zhenyu Wang <> at Jun 7th, 2021
#   Fix from python2 -> python3
# modified by Xiaoshan Luo <luoxs@calypso.cn> at Jun 8th, 2021
#   Add comment and change format style

import glob
import os

import numpy as np

# from sympy import EX


def pressure(outcar_path, contcar_path):
    info = os.popen(f"grep external {outcar_path} | tail -n 1").read().split()
    infokb = os.popen(f"grep 'in kB' {outcar_path} | tail -n 1").read().split()
    try:
        external = float(info[3])
        pullay = float(info[8])
        ex_x = float(infokb[2]) - float(pullay)
        ex_y = float(infokb[3]) - float(pullay)
        ex_z = float(infokb[4]) - float(pullay)
    except Exception:
        return False
    f = open(contcar_path)
    for i in range(0, 6):
        line = f.readline()
    try:
        natom = sum(map(int, line.split()))
    except Exception:
        line = f.readline()
        natom = sum(map(int, line.split()))
    f.close()
    f = open(outcar_path)
    numi = []
    fff = []
    for line in f:
        fff.append(line)
    for i in range(len(fff)):
        if 'TOTAL-FORCE' in fff[i]:
            numi.append(i)
    f.close()
    atomf = []
    latom = True
    for i in range(int(natom)):
        atomf.append(list(map(float, fff[i + numi[-1] + int(2)].split())))
    for i in range(int(natom)):
        if abs(atomf[i][3]) < 0.3 and abs(atomf[i][4]) < 0.3 and abs(atomf[i][5]) < 0.3:
            latom = True
        else:
            latom = False
            return False
    lstress = True
    if abs(pullay - 0.0) < 1.0:
        if (
            abs(ex_x) < 60.0
            and abs(ex_y) < 60.0
            and abs(ex_z) < 60.0
            and abs(external) < 20.0
        ):
            return True
        else:
            return False
    elif pullay < 1000.0:
        if (
            abs(ex_x) < 150.0
            and abs(ex_y) < 150.0
            and abs(ex_z) < 150.0
            and abs(external) < 50.0
        ):
            return True
        else:
            return False
    elif pullay < 10000.0:
        if (
            abs(ex_x) < 300.0
            and abs(ex_y) < 300.0
            and abs(ex_z) < 300.0
            and abs(external) < 100.0
        ):
            return True
        else:
            return False
    else:
        if (
            abs(ex_x) < 1500.0
            and abs(ex_y) < 1500.0
            and abs(ex_z) < 1500.0
            and abs(external) < 500.0
        ):
            return True
        else:
            return False


def outEnergy(outcar_path, contcar_path):
    f = open(contcar_path)
    for i in range(0, 6):
        line = f.readline()
    try:
        natom = sum(map(int, line.split()))
    except Exception:
        line = f.readline()
        natom = sum(map(int, line.split()))
    f.close()
    osawk = (
        '''awk '/enthalpy is  TOTEN    =/ {print $5}' '''
        + str(outcar_path.absolute())
        + ''' | tail -1'''
    )
    b = os.popen(osawk).read()
    a = b[0:-1]
    try:
        e = float(a) / natom
    except Exception:
        e = np.nan
    return e


def hardness():
    try:
        f = open('hardnes.dat')
    except FileNotFoundError:
        return 0
    try:
        b = f.readline()
    except Exception:
        return 0
    try:
        h = float(b)
    except Exception:
        return 0
    return h


def get_enthalpy(outcar_path, contcar_path):
    if pressure(outcar_path, contcar_path):
        e = outEnergy(outcar_path, contcar_path)
    else:
        # e = np.nan
        e = 610612509
    # h = hardness()
    return e


if __name__ == '__main__':
    pass
    # print(enthalpy('./'))
