from collections import defaultdict
from typing import List, Sequence
from pprint import pprint

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

def find_cless_ops(cless:List[np.ndarray], spgops:List[SymmOp])-> List[np.ndarray]:
    """
    在一组对称操作中查找与给定类中的矩阵相匹配的操作。

    Args:
        cless: 给定类，包含多个矩阵的列表。
        spgops: 包含多个对称操作的列表，每个操作都是一个矩阵。

    Returns:
        匹配的对称操作列表opslist: List[SymmOp]。
    """
    opslist = []
    for rot in cless:
        for op in spgops:
            if np.allclose(op.rotation_matrix, rot):
                opslist.append(op)
                break
    return opslist

def get_rot_name(rot:np.ndarray)-> str: ### 这个函数废除了，不用了!!!, get_operation_type代替了它的功能
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

def get_operation_type(Lattice, rotation):
    """
    Calculates the rotation axis and angle of the symmetry and if it 
    preserves handedness or not. 计算对称的旋转轴和角度，以及它是否保留手性。

    Returns
    -------
    tuple
        The first element is an array describing the rotation axis. The 
        second element describes the rotation angle. The third element is a 
        boolean, `True` if the symmetry preserves handedness 
        (determinant -1).
    """
    rotxyz = Lattice.T.dot(rotation).dot(np.linalg.inv(Lattice).T) #这里是为什么？？？？似乎是将旋转矩阵从晶格坐标转化为直角坐标
    # print ("rotation in real space:\n",rotxyz)
    E, V = np.linalg.eig(rotxyz)
    if not np.isclose(abs(E), 1).all():
        raise RuntimeError(
            "some eigenvalues of the rotation are not unitary")
            # 旋转的一些特征值是不酉的, 即, S·S-1 != E
    if E.prod() < 0:
        inversion = True
        E *= -1
    else:
        inversion = False
    idx = np.argsort(E.real)
    E = E[idx]
    V = V[:, idx]
    axis = V[:, 2].real
    if np.isclose(E[:2], 1).all():
        angle = 0
    elif np.isclose(E[:2], -1).all():
        angle = np.pi
    else:
        angle = np.angle(E[0])
        v = V[:, 0]
        s = np.real(np.linalg.det([v, v.conj(), axis]) / 1.j)
        if np.isclose(s, -1):
            angle = 2 * np.pi - angle
        elif not np.isclose(s, 1):
            raise RuntimeError("the sign of rotation should be +-1")
    return (axis, angle, inversion)

def get_cless_name(cless:List[SymmOp], lattice: np.ndarray)-> str:
    """
    根据一组共轭的对称操作推断类的共轭类名称。

    Args:
        cless: 一个共轭类列表包含多个对称操作的列表。

    Returns:
        共轭类的类型名称。
    """
    axis, angle, inversion = get_operation_type(struct.lattice.matrix, cless[0].rotation_matrix)
    axis_string  = list(map(str, np.round(axis, decimals=2)))
    angle_degree = str(np.round(np.degrees(angle), decimals=0))
    if inversion:
        name = angle_degree + '_' + "inversion" + '_axis__' + '_'.join(axis_string)
    else:
        name = angle_degree + '_axis__' + '_'.join(axis_string)
    
    return name

def classify_ops(spgops:List[SymmOp], lattice: np.ndarray) -> dict[str:List[SymmOp]]:
    """
    对称操作按照共轭类分类函数。

    Args:
        spgops: 包含多个对称操作的列表。

    Returns:
        按照共轭类分类后的对称操作字典，键是类的名称，值是包含多个对称操作的列表。
    """
    clesseslist: List[array] = []
    for idx, op in enumerate(spgops):
        # print(op.as_xyz_str())
        bm = [np.dot(np.dot(np.linalg.inv(x.rotation_matrix), op.rotation_matrix), x.rotation_matrix) for x in spgops]
        # print(bm)
        cless = kill_duplicated_element(bm) # 删除重复的矩阵，得到一个共轭类
        # print(cless)
        # for m in cless:
        #     print(detailed_info(m))
        # break
        clesseslist.append(cless) # 将类添加到类列表中

    clesseslist = kill_duplicated_clesses(clesseslist) # 删除重复的共轭类，得到一个共轭类
    # print((len(clesseslist)))
    clessesdict = defaultdict(list)
    for idx, cless in enumerate(clesseslist):
        ops:List[SymmOp] = find_cless_ops(cless, spgops=spgops)
        clessname:str = get_cless_name(ops, lattice)
        clessesdict[f'{clessname}'] = ops

    clessesdict:dict[List[SymmOp]] = dict(clessesdict)

    # if len(clesseslist) != len(clessesdict):
    #     raise ValueError('The size of clesseslist {} != the size of clessesdict {}'.format(len(clesseslist), len(clessesdict)))

    return clessesdict

def generate_mapping_matrix(mapping:dict[int:int], size:int)-> List[List[int]]:
    """
    通过给定的映射关系生成矩阵 A

    Args:
        mapping: 一个字典，表示映射关系，键为第一列的元素，值为第二列的元素

    Returns:
        矩阵 A
    """
    # 确定矩阵 A 的大小
    A = [[0] * size for _ in range(size)]
    # print(A)
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
    
    Lattice = struct.lattice.matrix
    rotation = op.rotation_matrix
    mapping = {}
    for idx1, site1 in enumerate(struct):
        
        # 方式一
        newcoords = op.operate(site1.frac_coords)
        newsite = PeriodicSite(species=site1.species, coords=newcoords, to_unit_cell=True, lattice=struct.lattice, coords_are_cartesian=False)
        # print(site1.frac_coords, newsite.frac_coords)
        # s = np.array(([
        #     [1, -np.sin(np.pi/6), 0],
        #     [0,  np.cos(np.pi/6), 0],
        #     [0,  0,  1],
        # ]))
        # s_1 = np.linalg.inv(s)
        # new_rotation_matrix = np.dot(np.dot(s_1, op.rotation_matrix), s)
        # new_rotation_matrix = np.dot(np.dot(s, op.rotation_matrix), s_1)
         # pprint(op.rotation_matrix)
        # pprint(new_rotation_matrix)
        
        # 方式二
        # rotxyz = Lattice.T.dot(rotation).dot(np.linalg.inv(Lattice).T) 
        # newcoords = np.dot(rotxyz, site1.coords)
        # newsite = PeriodicSite(species=site1.species, coords=newcoords, to_unit_cell=True, lattice=struct.lattice, coords_are_cartesian=True)
        # print(site1.frac_coords, newsite.frac_coords)

        for idx2, site2 in enumerate(struct):
            # print("判断 new-{} {} 是否与 {} {} 等价".format(idx1, newsite.frac_coords, idx2, site2.frac_coords))
            # print("{} {}-> {} {}".format(idx1, site1, idx2, newsite))
            # input()
            # print(np.allclose(newsite.frac_coords, site2.frac_coords, rtol=1e-3))
            if np.allclose(newsite.frac_coords, site2.frac_coords, rtol=1e-3):  # pymatgen 源代码检查了元素是否相等，坐标是否相等，性质是否相等 self.species == other.species and np.allclose(self.coords, other.coords, atol=Site.position_atol) and self.properties == other.properties
                mapping[idx1] = idx2
                break
        else:
            print("没有找到经过对称操作 {} 作用后，与{} {}等价的原子".format(op, idx1, site1.frac_coords))
        # print(mapping)
    # print(mapping)
    return mapping

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
    # print(ops)
    for op in ops:
        # print(op)
        mapping = get_atoms_mapping(op, struct)
        # pprint(mapping)
        A = generate_mapping_matrix(mapping, size=len(struct))
        trc = np.trace(A)
        trcs.append(trc)
        # pprint(A)
    if len(set(trcs)) == 1:
        return trc, A
    else:
        raise ValueError("The value of traces in the same conjugate class is different")

def get_atomic_character(new_spgops:dict[str:SymmOp], struct:Structure)-> tuple[dict[str:int], List[List[int]]]: 
    """
    获取原子的特征标

    Args:
        new_spgops: 包含所有分好共轭类的字典，键是共轭类的名称，值是相应的对称操作列表。
        struct: 结构。

    Returns:
        包含原子特征的字典，键是操作类的名称，值是对应的特征值。
    """
    atomic_charater = {}
    atomic_mapmatrix= {}
    for name, ops in new_spgops.items():
        # print(name)
        trc, matrix = get_conjugate_character(ops, struct)
        # print(trc, matrix)
        # input()
        atomic_charater[name] = trc
        atomic_mapmatrix[name] = matrix
    return (atomic_charater, atomic_mapmatrix)

if __name__ == '__main__':
    struct = Structure.from_file('../test/POSCAR_CaH6')
    partial_struct = Structure.from_file("../test/POSCAR_CaH6_partical")
    # spgops = spg.get_space_group_operations()
    spg = SpacegroupAnalyzer(struct)
    spgops = spg.get_point_group_operations()
    new_spgops = classify_ops(spgops, partial_struct.lattice)
    conjugate_charater, atomic_mapmatrix = get_atomic_character(new_spgops, partial_struct)

    # names = ['identity_1','rotation90_6', 'rotation180_3', 'rotation120_8', 'rotation180_6',  'inversion_1', 'rotoinversion90_6', 'rotoinversion180_3',  'rotoinversion120_8', 'rotoinversion180_6', ]
    names = [
        '0.0_axis__0.0_0.99_0.13',
        '90.0_axis__0.0_0.0_1.0',
        '180.0_axis__0.0_0.0_1.0',
        '120.0_axis__0.58_0.58_0.58',
        '180.0_axis__-0.71_0.71_0.0',
        '0.0_inversion_axis__0.0_0.99_-0.13',
        '450.0_inversion_axis__0.0_0.0_1.0',
        '180.0_inversion_axis__0.0_0.0_1.0',
        '-120.0_inversion_axis__-0.58_-0.58_-0.58',
        '180.0_inversion_axis__0.71_-0.71_-0.0',
                ]
    mult  = np.array([1,6,3,8,6, 1,6,3,8,6])
    red_reps = np.array([conjugate_charater[name] for name in names])


    irred_repss={
        'A1g' : np.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1]),
        'A1u' : np.array([1, 1, 1, 1, 1,-1,-1,-1,-1,-1]),
        'A2g' : np.array([1,-1, 1, 1,-1, 1,-1, 1, 1,-1]),
        'A2u' : np.array([1,-1, 1, 1,-1,-1, 1,-1,-1, 1]),
        'Eg'  : np.array([2, 0, 2,-1, 0, 2, 0, 2,-1, 0]),
        'Eu'  : np.array([2, 0, 2,-1, 0,-2, 0,-2, 1, 0]),
        'T2u' : np.array([3,-1,-1, 0, 1,-3, 1, 1, 0,-1]),
        'T2g' : np.array([3,-1,-1, 0, 1, 3,-1,-1, 0, 1]),
        'T1u' : np.array([3, 1,-1, 0,-1,-3,-1, 1, 0, 1]),
        'T1g' : np.array([3, 1,-1, 0,-1, 3, 1,-1, 0,-1]),
        }
    for irred_reps_name, irred_reps_character in irred_repss.items():
        x=np.sum(mult*irred_reps_character*(np.array(red_reps)))/np.sum(mult)
        if x != 0:
            print(irred_reps_name)