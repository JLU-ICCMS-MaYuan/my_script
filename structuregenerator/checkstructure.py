import logging
from itertools import product
from typing import *

import numpy as np

from ase.atoms import Atoms
from ase.io import read
from ase.build import make_supercell

from sklearn.preprocessing import minmax_scale

logging.basicConfig(level = logging.INFO,format = '%(asctime)s|%(name)s|%(levelname)s|%(message)s')
logger = logging.getLogger(__name__)

def MinMax_normalize(data: np.ndarray):
    '''
    离差标准化
    x = x - min / max - min
    '''
    max = data.max()
    min = data.min()
    result = (data - min) / (max -min)
    return result

def MeanVariance_normalize(X: np.ndarray):
    '''
    计算均值方差归一化
    把所有数据归一到均值为0方差为1的数据中
    既适用于数据没有明显的边界，有可能存在极端数据值的情况。也适用于数据有明显边界的情况
    x_scale =  ( x - x_mean ) / s
        s      是一维数组x的标准差
        x_mean 是一维数组的平均值
        x_scale是一维数组的均值方差归一化值
    '''

    result = (X[:] - np.mean(X[:])) / np.std(X[:])

    # np.std 计算标准差 
    #   平均值 = x_mean = np.average((x1 + x2 + ... + xn))
    #   标准差 = [(x1 - x_mean) + (x2 - x_mean) + ... + (xn - x_mean)] / n
    return result

def checkHDistance(struc: Atoms):
    for a in struc:
        for b in struc:
            if (a.index != b.index) and (a.symbol == "H") and (b.symbol == "H"):
                d = struc.get_distance(a.index, b.index)
                if d < 0.8:
                    return False
    else:
        return True

def checkHCage(struc: Atoms):
    lower_limit = 20 # 构成一个氢笼子需要氢原子的下限。只要达到了这个下限，就认为该非氢原子周围围满了足够多的氢原子用于构成氢笼子
    for a in struc:
        if (a.symbol != "H"):
            contribution_for_cage = 0
            for b in struc:
                if(b.symbol == "H"):
                    d = struc.get_distance(a.index, b.index)
                    if (1.8 <= d) and (d <= 2.2) :
                        contribution_for_cage += 1
                    if d < 1.8:
                        logger.info(f"non-H element {a.symbol}{a.index} is too closed with the hydrogen {b.symbol}{b.index}")
                        return False
            else:
                if contribution_for_cage >= lower_limit:
                    logger.info(f"There are    enough hydrogens aroud {a.symbol}{a.index}, whose amount is {contribution_for_cage}")
                else:
                    logger.info(f"There are no enough hydrogens aroud {a.symbol}{a.index}, whose amount is {contribution_for_cage}!")
    else:
        return False

def checkHdensity(struc: Atoms):
    # 1. 扩包
    # 2. 计算晶格实空间内氢原子的密度：
    #    把晶格打网格，计算每一个网格点在某平均内平均氢原子的个数，即在这一点的H的密度为
    #       Hdensity = V_ball / N_Hatoms
    # 3. 评估这些点所对应的氢的密度有没有特别大的点和特别小的点，使得他们整体差异非常大
    #       Hdensity(x, y, z) = value
    positions : np.ndarray = struc.positions
    xmin, ymin, zmin = positions.min(axis=0)
    xmax, ymax, zmax = positions.max(axis=0)

    x_mesh = np.arange(xmin+1.0, xmax, 1.0, dtype=np.float64)
    y_mesh = np.arange(ymin+1.0, ymax, 1.0, dtype=np.float64)
    z_mesh = np.arange(zmin+1.0, zmax, 1.0, dtype=np.float64)

    points = product(x_mesh, y_mesh, z_mesh)
    __coords_Hdensity = []
    __total_Hdensity = []
    R = 2.0
    for p in points:
        Nhydrogen = 0
        p_Hdensity = 0
        for a in struc:
            if a.symbol == "H":
                d = np.linalg.norm(p-a.position)
                if d < R:
                    Nhydrogen += 1

        if not np.isclose(Nhydrogen, 0.0):
            p_Hdensity = (4/3 * np.pi * np.power(R, 3)) / Nhydrogen
            __coords_Hdensity.append(np.concatenate((p, [p_Hdensity]), axis=0))
            __total_Hdensity.append(p_Hdensity)

    __total_Hdensity = np.array(__total_Hdensity)[:, np.newaxis]
    __coords_Hdensity = np.array(__coords_Hdensity)

    _total_Hdensity = minmax_scale(
        X=__total_Hdensity,      # 输入一个矩阵，即一个二维数组，如果你输入一维数组，会报错的, 你需要自行扩维
        feature_range=(0, 1),   # 设置将数据归一化到哪个范围内
        axis=0,                 # 以每一列为一个整体，计算每一列的数据对应的标准差标准化的值
        copy=True,              # 这个开关：布尔值，是否拷贝一份数据以避免在原数据上进行操作，默认为True
    )

    _coords_Hdensity = minmax_scale(
        X=__coords_Hdensity[:, 3], # 输入一个矩阵，即一个二维数组，如果你输入一维数组，会报错的, 你需要自行扩维
        feature_range=(0, 1),   # 设置将数据归一化到哪个范围内
        axis=0,                 # 以每一列为一个整体，计算每一列的数据对应的标准差标准化的值
        copy=True,              # 这个开关：布尔值，是否拷贝一份数据以避免在原数据上进行操作，默认为True
    )

    index_of_not1and0 = np.where((~np.isclose(_total_Hdensity, 1.0)) & (~np.isclose(_total_Hdensity, 0.0)))
    total_Hdensity = _total_Hdensity[index_of_not1and0]
    std_of_totalHdensity = np.std(total_Hdensity ,axis=0)
    logger.info(std_of_totalHdensity)
    return np.array(total_Hdensity) 



if __name__ == "__main__":
    # struc = read(f"/Users/macbookpro/my_code/my_script/test/LaH10/LaH10.vasp")
    struc = read(f"/Users/macbookpro/Library/CloudStorage/OneDrive-mails.jlu.edu.cn/氢化物结构/CaH6.vasp")
    dim = [[2, 0, 0],[0, 2, 0],[0, 0, 2]]
    superstruc = make_supercell(struc, dim)
    print(len(superstruc))
    if checkHDistance(superstruc):
        print("step 1st checkHDistance succeeded !")
        if checkHCage(superstruc):
            print("step 2nd checkHCage succeeded !")
            coords_Hdensity = checkHdensity(superstruc) 
        else:
            print("step 2nd checkHCage failed !")
    else:
        print("step 1st checkHDistance failed !")




