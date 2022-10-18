import logging
from itertools import product
from collections import defaultdict
from typing import *

import numpy as np

from ase.atoms import Atoms
from ase.io import read
from ase.build import make_supercell

from pymatgen.core.structure import Structure
from pymatgen.io.ase import AseAtomsAdaptor

from sklearn.preprocessing import minmax_scale

logging.basicConfig(level = logging.INFO,format = '%(asctime)s|%(name)s|%(levelname)s|%(message)s')
logger = logging.getLogger(__name__)


def checkHDistance(pmg_struc, center_indices, points_indices, images, distances):
    for center, points, imgs, dist in zip(center_indices, points_indices, images, distances):
        if (str(pmg_struc.sites[center].specie) == "H") and \
           (str(pmg_struc.sites[points].specie) == "H"):
            if dist < 0.8:
                return False
    else:
        return True

def checkHCage(pmg_struc, center_indices, points_indices, images, distances):

    large_size = defaultdict(int)
    for center, points, imgs, dist in zip(center_indices, points_indices, images, distances):
        if str(pmg_struc.sites[center].specie) != "H" and 1.7 <= dist <= 2.2 :
            large_size[str(pmg_struc.sites[center].specie)] += 1

    small_size = defaultdict(int)
    for center, points, imgs, dist in zip(center_indices, points_indices, images, distances):
        if str(pmg_struc.sites[center].specie) != "H" and 1.3 <= dist <= 1.8 :
            small_size[str(pmg_struc.sites[center].specie)] += 1

    el_amt = pmg_struc.composition.get_el_amt_dict()
    specie_without_H = list(el_amt.keys())
    specie_without_H.remove("H")
    cage_size_without_reducing = {}
    for key in specie_without_H:
        from_big_cage = large_size.get(key, 0)
        from_small_cage = small_size.get(key, 0)
        if from_big_cage > from_small_cage:
            cage_size_without_reducing[key] = from_big_cage
        else:
            cage_size_without_reducing[key] = from_small_cage

    cage_size_with_reducing = defaultdict(int)
    for key, value in cage_size_without_reducing.items():
        if key in list(el_amt.keys()):
            around = value / el_amt[key]
            cage_size_with_reducing[key] = around

    cage_size_with_reducing = dict(cage_size_with_reducing)
    for key, value in cage_size_with_reducing.items():
        if value < 8:
            return False
    else:
        return cage_size_with_reducing
    
def checkHdensity(pmg_struc, center_indices, points_indices, images, distances):
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
        X=__total_Hdensity,     # 输入一个矩阵，即一个二维数组，如果你输入一维数组，会报错的, 你需要自行扩维
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

def check(pmg_struc: Structure):

    (
        center_indices,
        points_indices,
        images,
        distances,
    ) = pmg_struc.get_neighbor_list(
        r=2.2, 
        sites=pmg_struc.sites, 
        numerical_tol=1e-8
    )
    if checkHDistance(pmg_struc, center_indices, points_indices, images, distances):
        res: Dict = checkHCage(pmg_struc, center_indices, points_indices, images, distances)
        if res:
            logger.info(f"finally successed! The program find some cages which is {res}")
            return True
        else:
            # logger.info("check H-cage failed")
            return False
    else:
        # logger.info("check H-distance failed")
        return False
    

if __name__ == "__main__":

    from pathlib import Path
    # struc = read(f"/Users/macbookpro/my_code/my_script/test/LaH10/LaH10.vasp")
    # struc = read(f"/Users/macbookpro/Library/CloudStorage/OneDrive-mails.jlu.edu.cn/氢化物结构/CaH6.vasp")
    # struc = read("/home/mayuan/mycode/my_script/test/clathrate/CaH6.vasp")
    # struc = read("/home/mayuan/mycode/my_script/test/clathrate/LaH10.vasp")
    # struc = read("/home/mayuan/mycode/my_script/test/clathrate/LaBH8.vasp")
    # struc = read("/home/mayuan/mycode/my_script/test/clathrate/Li2MgH16.vasp")
    # struc = read("/home/mayuan/mycode/my_script/test/clathrate/ScH9.vasp")
    # struc = read("/home/mayuan/mycode/my_script/test/clathrate/CaYH12.vasp")
    
    # for path in Path("/home/mayuan/mycode/my_script/test/194/unstable_structs").glob("UCell_*"):
    for path in Path("/work/home/may/calypso-prediction/194/result/unstable_structs").glob("UCell_*"):
        struc = read(path)
        res = check(struc)
        if res:
            print(struc.symbols, path)




