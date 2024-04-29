from collections import defaultdict
from typing import List

import numpy as np
from numpy import array

from pymatgen.symmetry.analyzer import SymmOp, SpacegroupAnalyzer
from pymatgen.core.structure import Structure

def is_in_matrix(m:array, ms:List[array]) -> bool:
    """
    判断矩阵 m 是否在矩阵组成的列表 ms 中。

    Args:
        m: 要判断的矩阵。
        ms: 包含多个矩阵的列表。

    Returns:
        如果矩阵 m 存在于列表 ms 中，则返回 True，否则返回 False。
    """
    flag = True  # 初始化标志位为 True
    for _ in ms: # 遍历列表 ms 中的每个矩阵
        if np.allclose(m, _):  # 使用 np.allclose() 判断两个矩阵是否相似
            flag = True  # 如果找到相似的矩阵，则将标志位设为 True
            break  # 找到相似矩阵后退出循环
    else:  # 如果未找到相似的矩阵
        flag = False  # 将标志位设为 False
    return flag  # 返回标志位

def kill_duplicated_element(ms: List[np.ndarray]) -> List[np.ndarray]:
    """
    在一组 3*3 的矩阵中，找到重复的并删除。

    Args:
        ms: 包含多个 3*3 的矩阵的列表。

    Returns:
        删除重复矩阵后的新列表。
    """
    new_ms = []  # 初始化新列表
    for m in ms:  # 遍历原始列表中的每个矩阵
        if not new_ms:  # 如果新列表为空
            new_ms.append(m)  # 直接将当前矩阵添加到新列表中
        else:
            if not is_in_matrix(m=m, ms=new_ms):  # 如果当前矩阵不在新列表中
                new_ms.append(m)  # 将当前矩阵添加到新列表中
    return new_ms  # 返回新列表

def is_equal_for_cless(cl1: List[np.ndarray], cl2: List[np.ndarray]) -> bool:
    """
    判断两个类是否相同。

    Args:
        cl1: 第一个类，包含多个矩阵的列表。
        cl2: 第二个类，包含多个矩阵的列表。

    Returns:
        如果两个类相同（即包含相同的矩阵），返回 True，否则返回 False。
    """
    flag = True  # 初始化标志位为 True
    for m1 in cl1:  # 遍历第一个类中的每个矩阵
        if is_in_matrix(m1, cl2):  # 调用 is_in_matrix() 函数判断矩阵是否在第二个类中
            flag = True  # 如果找到相同的矩阵，则将标志位设为 True
            break  # 找到相同的矩阵后退出循环
    else:  # 如果第一个类中的所有矩阵都不在第二个类中
        flag = False  # 将标志位设为 False
    return flag  # 返回标志位

def is_in_clesses(cl1: List[np.ndarray], clesses: List[List[np.ndarray]]) -> bool:
    """
    判断一个类是否在一组类中。

    Args:
        cl1: 要判断的类，包含多个矩阵的列表。
        clesses: 包含多个类的列表，每个类都是包含多个矩阵的列表。

    Returns:
        如果要判断的类在一组类中，则返回 True，否则返回 False。
    """
    flag = True  # 初始化标志位为 True
    for cl2 in clesses:  # 遍历一组类中的每个类
        if is_equal_for_cless(cl1, cl2):  # 调用 is_equal_for_cless() 函数判断两个类是否相同
            flag = True  # 如果找到相同的类，则将标志位设为 True
            break  # 找到相同的类后退出循环
    else:  # 如果要判断的类不在一组类中
        flag = False  # 将标志位设为 False
    return flag  # 返回标志位

def kill_duplicated_clesses(clesses:List[List[np.ndarray]]) -> List[List[np.ndarray]]:
    """
    删除一组类中的重复类。

    Args:
        clesses: 包含多个类的列表，每个类都是包含多个矩阵的列表。

    Returns:
        删除重复类后的新列表。
    """
    new_clesses = []
    for cl in clesses: # cl此时是一个类，里面有很多个3*3的矩阵
        if not new_clesses:
            new_clesses.append(cl)
        else:
            if not is_in_clesses(cl1=cl, clesses=new_clesses):
                new_clesses.append(cl)
    return new_clesses

def find_cless_ops(cless:List[np.ndarray])-> List[np.ndarray]:
    """
    在一组对称操作中查找与给定类中的矩阵相匹配的操作。

    Args:
        cless: 给定类，包含多个矩阵的列表。
        pgops: 包含多个对称操作的列表，每个操作都是一个矩阵。

    Returns:
        匹配的对称操作列表tmplist: List[SymmOp]。
    """
    tmplist = []
    for rot in cless:
        for op in pgops:
            if np.allclose(op.rotation_matrix, rot):
                tmplist.append(op)
                break
    return tmplist

def get_rot_name(rot:np.ndarray)-> str:
    """
    根据旋转矩阵推断对称操作的类型。

    Args:
        rot: 旋转矩阵。

    Returns:
        对称操作的类型。
    """
    typename = ''
    det = np.linalg.det(rot) 
    if det > 0: # det(rot)=1: rotation; 
        trc = np.trace(rot)
        theta = np.degrees(np.arccos((trc-1)/2))
        # print(trc); print(theta)
        if np.allclose(theta, 0.0):
            typename = "identity"
        else:
            typename = "rotation"+str(int(theta))
    else: # det(rot)=-1: rotoinversion
        trc = -np.trace(rot)
        theta = np.degrees(np.arccos((trc-1)/2))
        # print(trc); print(theta)
        if np.allclose(theta, 0.0):
            typename = "inversion"
        else:
            typename = "rotoinversion"+str(int(theta))

    return typename

def get_cless_name(cless:List[SymmOp])-> str:
    """
    根据一组共轭的对称操作推断类的共轭类名称。

    Args:
        cless: 一个共轭类列表包含多个对称操作的列表。

    Returns:
        共轭类的类型名称。
    """
    namelist = []
    for op in cless:
        name = get_rot_name(op.rotation_matrix)
        namelist.append(name)
    
    if len(set(namelist)) == 1:
        return namelist[0]+ "_" + str(len(namelist))
    else:
        raise ValueError("Error: namelist={}".format(namelist))

def classify_ops(pgops:List[SymmOp]) -> dict[str:List[SymmOp]]:
    """
    对称操作按照共轭类分类函数。

    Args:
        pgops: 包含多个对称操作的列表。

    Returns:
        按照共轭类分类后的对称操作字典，键是类的名称，值是包含多个对称操作的列表。
    """
    clesseslist: List[array] = []
    for idx, op in enumerate(pgops):
        # print(op.as_xyz_str())
        bm = [np.dot(np.dot(np.linalg.inv(x.rotation_matrix), op.rotation_matrix), x.rotation_matrix) for x in pgops]
        # print(bm)
        cless = kill_duplicated_element(bm) # 删除重复的矩阵，得到一个共轭类
        # print(cless)
        # for m in cless:
        #     print(detailed_info(m))
        # break
        clesseslist.append(cless) # 将类添加到类列表中

    clesseslist = kill_duplicated_clesses(clesseslist) # 删除重复的共轭类，得到一个共轭类

    clessesdict = defaultdict(list)
    for idx, cless in enumerate(clesseslist):
        tmp:List[SymmOp] = find_cless_ops(cless)
        clessname:str = get_cless_name(tmp)
        clessesdict[f'{clessname}'] = tmp

    clessesdict:dict[List[SymmOp]] = dict(clessesdict)
    return clessesdict

def generate_mapping_matrix(mapping)-> List[List[int]]:
    """
    通过给定的映射关系生成矩阵 A

    Args:
        mapping: 一个字典，表示映射关系，键为第一列的元素，值为第二列的元素

    Returns:
        矩阵 A
    """
    # 确定矩阵 A 的大小
    size = len(mapping)
    A = [[0] * size for _ in range(size)]

    # 遍历映射关系，填充矩阵 A
    for idx, key in enumerate(mapping):
        value = mapping[key]
        A[idx][value] = 1
    # print(mapping)
    return A


def get_atoms_mapping(op:SymmOp, struct:Structure)-> dict[int:int]:
    """
    根据对称操作将原子映射到新的原子位置。

    Args:
        op: 对称操作。
        struct: 原始晶体结构。

    Returns:
        映射字典，键是原始原子的索引，值是新原子的索引。
    """
    from pymatgen.core.sites import PeriodicSite
    
    mapping = {}
    for idx1, site1 in enumerate(struct):
        newcoords = op.operate(site1.frac_coords)
        newsite = PeriodicSite(species=site1.species, coords=newcoords, to_unit_cell=True, lattice=struct.lattice)
        for idx2, site2 in enumerate(struct):
            if np.allclose(newsite.coords, site2.coords, rtol=1e-3):  # pymatgen 源代码检查了元素是否相等，坐标是否相等，性质是否相等 self.species == other.species and np.allclose(self.coords, other.coords, atol=Site.position_atol) and self.properties == other.properties
                # print("{} {}-> {} {}".format(idx1, site1.frac_coords, idx2, newsite.frac_coords))
                # print("{} {}-> {} {}".format(idx1, site1, idx2, newsite))
                mapping[idx1] = idx2
    return mapping

# for name in ['identity_1', 'rotation180_3', 'rotation180_6', 'rotation120_8', 'rotation90_6', 'inversion_1', 'rotoinversion180_3', 'rotoinversion180_6', 'rotoinversion120_8', 'rotoinversion90_6']:
#     print(get_atoms_mapping(new_pgops[name][0], struct=struct))

def get_conjugate_character(ops:List[SymmOp], struct:Structure)-> int:
    """
    获取给定共轭类的特征标

    Args:
        ops: 共轭类，包含多个对称操作的列表。
        struct: 结构。

    Returns:
        共轭类的特征标
    """
    trcs = []
    for op in ops:
        mapping = get_atoms_mapping(op, struct)
        A = generate_mapping_matrix(mapping)
        trc = np.trace(A)
        trcs.append(trc)

    if len(set(trcs)) == 1:
        return trc, A

# for name in ['identity_1', 'rotation180_3', 'rotation180_6', 'rotation120_8', 'rotation90_6', 'inversion_1', 'rotoinversion180_3', 'rotoinversion180_6', 'rotoinversion120_8', 'rotoinversion90_6']:
#     print(get_conjugate_character(new_pgops[name], struct=struct))

def get_atomic_character(new_pgops:dict[str:SymmOp], struct:Structure)-> dict[str:int]: 
    """
    获取原子的特征标

    Args:
        new_pgops: 包含所有分好共轭类的字典，键是共轭类的名称，值是相应的对称操作列表。
        struct: 结构。

    Returns:
        包含原子特征的字典，键是操作类的名称，值是对应的特征值。
    """
    atomic_charater = {}
    atomic_mapmatrix= {}
    for name, ops in new_pgops.items():
        trc, matrix = get_conjugate_character(ops, struct)
        atomic_charater[name] = trc
        atomic_mapmatrix[name] = matrix
    return atomic_charater, atomic_mapmatrix



if __name__ == "__main__":

    from pprint import pprint

    struct = Structure.from_file("../test/POSCAR_LaH10")
    spg = SpacegroupAnalyzer(struct)
    pgops = spg.get_point_group_operations()
    new_pgops = classify_ops(pgops)
    conjugate_charater, atomic_mapmatrix = get_atomic_character(new_pgops, struct)
    # print(conjugate_charater)

    names = ['identity_1','rotation90_6', 'rotation180_3', 'rotation120_8', 'rotation180_6',  'inversion_1', 'rotoinversion90_6', 'rotoinversion180_3',  'rotoinversion120_8', 'rotoinversion180_6']
    for name in names:
        print(name)
        pprint(atomic_mapmatrix[name])
        pprint(conjugate_charater[name])
    red_reps = np.array([conjugate_charater[name] for name in names])-1
    print(red_reps)
    mult  = np.array([int(name.split('_')[-1]) for name in names])
    # print(mult)

    irred_repss={
        'A1g' : [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
        'A1u' : [1, 1, 1, 1, 1,-1,-1,-1,-1,-1],
        'A2g' : [1,-1, 1, 1,-1, 1,-1, 1, 1,-1],
        'A2u' : [1,-1, 1, 1,-1,-1, 1,-1,-1, 1],
        'Eg'  : [2, 0, 2,-1, 0, 2, 0, 2,-1, 0],
        'Eu'  : [2, 0, 2,-1, 0,-2, 0,-2, 1, 0],
        'T2u' : [3,-1,-1, 0, 1,-3, 1, 1, 0,-1],
        'T2g' : [3,-1,-1, 0, 1, 3,-1,-1, 0, 1],
        'T1u' : [3, 1,-1, 0,-1,-3,-1, 1, 0, 1],
        'T1g' : [3, 1,-1, 0,-1, 3, 1,-1, 0,-1],
        }

    for irred_reps_name, irred_reps_character in irred_repss.items():
        x=np.sum(np.array(mult)*np.array(irred_reps_character)*(np.array(red_reps)))/np.sum(mult)
        if x != 0:
            print(irred_reps_name)
    # from pymatgen.core.structure import Molecule
    # from pymatgen.symmetry.analyzer import PointGroupAnalyzer

    # water_coords = [(0.0, 0.0, 0.0), (0.757, 0.586, 0.0), (-0.757, 0.586, 0.0)]
    # water_species = ["H", "H", "O"]

    # # 创建水分子对象
    # water_molecule = Molecule(water_species, water_coords)
    # pg = PointGroupAnalyzer(water_molecule, matrix_tolerance=0.001)
    # pgops = pg.get_symmetry_operations()
    # print(pgops)
    # print(pg.get_pointgroup())
    # new_pgops:dict[List[SymmOp]] = classify_ops(pgops)
    # print(new_pgops.keys())

    # conjugate_charater = get_atomic_character(new_pgops, water_molecule)