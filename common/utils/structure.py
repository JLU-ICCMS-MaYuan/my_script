"""
结构文件处理工具

提供统一的结构文件读取、转换和操作功能。
支持POSCAR、CIF、XSF等格式。

作者：Claude
创建时间：2025-11-20
"""

from pathlib import Path
from typing import Dict, List, Tuple, Optional
import numpy as np


def read_poscar(poscar_file: Path) -> Dict:
    """
    读取VASP POSCAR文件

    Parameters
    ----------
    poscar_file : Path
        POSCAR文件路径

    Returns
    -------
    Dict
        包含结构信息的字典：
        - comment: 注释行
        - scaling: 缩放因子
        - lattice: 晶格矢量 (3x3)
        - species: 元素符号列表
        - numbers: 各元素原子数
        - coord_type: 坐标类型 ('Direct' or 'Cartesian')
        - positions: 原子位置列表
    """
    with open(poscar_file, 'r') as f:
        lines = f.readlines()

    # 解析注释
    comment = lines[0].strip()

    # 缩放因子
    scaling = float(lines[1].strip())

    # 晶格矢量
    lattice = []
    for i in range(2, 5):
        lattice.append([float(x) for x in lines[i].split()])
    lattice = np.array(lattice)

    # 元素符号（VASP 5.x格式）
    species_line = lines[5].strip().split()
    if species_line[0].isalpha():
        species = species_line
        numbers = [int(x) for x in lines[6].split()]
        coord_line_idx = 7
    else:
        # 旧格式，无元素符号
        species = None
        numbers = [int(x) for x in lines[5].split()]
        coord_line_idx = 6

    # 坐标类型
    coord_type = lines[coord_line_idx].strip()[0].upper()
    coord_type_str = 'Direct' if coord_type in ['D', 'd'] else 'Cartesian'

    # 原子位置
    positions = []
    total_atoms = sum(numbers)
    for i in range(coord_line_idx + 1, coord_line_idx + 1 + total_atoms):
        pos = [float(x) for x in lines[i].split()[:3]]
        positions.append(pos)

    return {
        'comment': comment,
        'scaling': scaling,
        'lattice': lattice,
        'species': species,
        'numbers': numbers,
        'coord_type': coord_type_str,
        'positions': positions,
    }


def write_poscar(structure: Dict, output_file: Path):
    """
    写入VASP POSCAR文件

    Parameters
    ----------
    structure : Dict
        结构信息字典（格式同read_poscar返回值）
    output_file : Path
        输出文件路径
    """
    with open(output_file, 'w') as f:
        # 注释
        f.write(f"{structure.get('comment', 'Structure')}\n")

        # 缩放因子
        f.write(f"{structure.get('scaling', 1.0):.10f}\n")

        # 晶格矢量
        lattice = structure['lattice']
        for row in lattice:
            f.write(f"  {row[0]:15.10f} {row[1]:15.10f} {row[2]:15.10f}\n")

        # 元素符号
        if structure.get('species'):
            f.write("  " + "  ".join(structure['species']) + "\n")

        # 各元素原子数
        f.write("  " + "  ".join(map(str, structure['numbers'])) + "\n")

        # 坐标类型
        f.write(f"{structure.get('coord_type', 'Direct')}\n")

        # 原子位置
        for pos in structure['positions']:
            f.write(f"  {pos[0]:15.10f} {pos[1]:15.10f} {pos[2]:15.10f}\n")


def poscar_to_qe_format(poscar_file: Path) -> Dict:
    """
    将POSCAR格式转换为QE pw.in格式

    Parameters
    ----------
    poscar_file : Path
        POSCAR文件路径

    Returns
    -------
    Dict
        QE格式的结构信息
    """
    structure = read_poscar(poscar_file)

    # QE格式
    qe_structure = {
        'ibrav': 0,  # 使用CELL_PARAMETERS
        'nat': sum(structure['numbers']),
        'ntyp': len(structure['species']) if structure['species'] else len(structure['numbers']),
        'lattice': structure['lattice'] * structure['scaling'],  # 应用缩放
        'species': structure['species'],
        'positions': structure['positions'],
    }

    return qe_structure


def get_lattice_parameters(lattice: np.ndarray) -> Tuple[float, float, float, float, float, float]:
    """
    从晶格矢量计算晶格参数

    Parameters
    ----------
    lattice : np.ndarray
        晶格矢量矩阵 (3x3)

    Returns
    -------
    Tuple[float, float, float, float, float, float]
        (a, b, c, alpha, beta, gamma)
        长度单位：与输入一致
        角度单位：度
    """
    a = np.linalg.norm(lattice[0])
    b = np.linalg.norm(lattice[1])
    c = np.linalg.norm(lattice[2])

    alpha = np.arccos(np.dot(lattice[1], lattice[2]) / (b * c))
    beta = np.arccos(np.dot(lattice[0], lattice[2]) / (a * c))
    gamma = np.arccos(np.dot(lattice[0], lattice[1]) / (a * b))

    # 转换为度
    alpha = np.degrees(alpha)
    beta = np.degrees(beta)
    gamma = np.degrees(gamma)

    return a, b, c, alpha, beta, gamma


def direct_to_cartesian(direct_coords: np.ndarray, lattice: np.ndarray) -> np.ndarray:
    """
    分数坐标转笛卡尔坐标

    Parameters
    ----------
    direct_coords : np.ndarray
        分数坐标 (N, 3)
    lattice : np.ndarray
        晶格矢量 (3, 3)

    Returns
    -------
    np.ndarray
        笛卡尔坐标 (N, 3)
    """
    return np.dot(direct_coords, lattice)


def cartesian_to_direct(cart_coords: np.ndarray, lattice: np.ndarray) -> np.ndarray:
    """
    笛卡尔坐标转分数坐标

    Parameters
    ----------
    cart_coords : np.ndarray
        笛卡尔坐标 (N, 3)
    lattice : np.ndarray
        晶格矢量 (3, 3)

    Returns
    -------
    np.ndarray
        分数坐标 (N, 3)
    """
    return np.dot(cart_coords, np.linalg.inv(lattice))
