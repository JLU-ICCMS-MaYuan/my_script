#!/usr/bin/env python3

import sys
import argparse
import numpy as np
from ase.io import read, write
from ase.build import make_supercell
from ase.build.tools import sort

# 自定义排序函数，按照指定顺序排序元素
def sort_by_custom_order(atoms, preferred_order):
    """
    Sort atoms according to a custom preferred order of elements.

    Parameters:
    atoms: ASE Atoms object
        The atomic structure to be sorted.
    preferred_order: list of str
        The desired order of elements, e.g., ['O', 'C', 'H'].

    Returns:
    ASE Atoms object
        The sorted atomic structure.
    """
    symbols = atoms.get_chemical_symbols()
    tags = [preferred_order.index(symbol) if symbol in preferred_order else len(preferred_order) for symbol in symbols]
    return sort(atoms, tags)

# 扩胞并按元素顺序排序
def expand_cell(input_structure, expansion_matrix, preferred_order=None):
    atoms = read(input_structure)
    
    # 生成扩胞结构
    supercell = make_supercell(atoms, expansion_matrix)
    
    # 按照用户指定的顺序进行排序（如果提供了顺序）
    if preferred_order:
        sorted_supercell = sort_by_custom_order(supercell, preferred_order)
    else:
        sorted_supercell = sort(supercell)  # 默认排序
    
    return sorted_supercell

# 解析命令行参数
def parse_arguments():
    parser = argparse.ArgumentParser(description='Generate supercells with optional diagonal or non-diagonal supercell matrices.')
    
    # 输入文件
    parser.add_argument('-i', '--input', type=str, default='POSCAR',
                        help='Input structure file (default: POSCAR)')

    # 对角扩胞设置
    parser.add_argument('-ds', '--diagonal_supercell', nargs=3, type=int,
                        help='Diagonal supercell expansion factors, e.g., 2 2 2 for a 2x2x2 diagonal supercell')

    # 非对角扩胞设置
    parser.add_argument('-nds', '--nondiagonal_supercell', nargs=9, type=int,
                        help='Non-diagonal supercell expansion matrix, e.g., 2 1 0 0 2 1 0 0 3 for a 3x3 non-diagonal supercell')

    # 元素排序顺序
    parser.add_argument('-o', '--order', nargs='+', type=str,
                        help='Custom element order for sorting, e.g., O C H for ordering oxygen first, followed by carbon and hydrogen')

    return parser.parse_args()

if __name__ == "__main__":
    info = '''python3 supercell_expansion.py -i POSCAR -ds 2 2 2
-ds 2 2 2 表示对角方向的 2x2x2 扩胞
python3 supercell_expansion.py -i POSCAR -nds 0 1 1 1 0 1 1 1 0
-nds 0 1 1 1 0 1 1 1 0 生成一个 3x3 的非对角扩胞
python3 supercell_expansion.py -i POSCAR -ds 2 2 2 -o O C H
-o O C H 表示按氧、碳、氢的顺序排序
'''

    # 解析命令行参数
    args = parse_arguments()

    # 读取输入文件
    input_structure = args.input

    # 检查是否提供了对角或非对角扩胞矩阵
    if args.diagonal_supercell:
        # 对角扩胞，生成对角扩胞矩阵
        supercell_matrix = np.diag(args.diagonal_supercell)
        print(f"Using diagonal supercell expansion matrix:\n{supercell_matrix}")
    elif args.nondiagonal_supercell:
        # 非对角扩胞，生成非对角扩胞矩阵
        supercell_matrix = np.array(args.nondiagonal_supercell).reshape((3, 3))
        print(f"Using non-diagonal supercell expansion matrix:\n{supercell_matrix}")
    else:
        print("Error: You must provide either a diagonal or non-diagonal supercell matrix.")
        sys.exit(1)

    # 检查是否提供了元素排序顺序
    if args.order:
        preferred_order = args.order
        print(f"Using custom element order: {preferred_order}")
    else:
        preferred_order = None

    # 生成扩胞并排序
    supercell = expand_cell(input_structure, supercell_matrix, preferred_order)

    # 输出扩胞后的结构文件
    output_file = 'POSCAR_supercell'
    write(output_file, supercell)
    print(f"Supercell structure written to {output_file}")
