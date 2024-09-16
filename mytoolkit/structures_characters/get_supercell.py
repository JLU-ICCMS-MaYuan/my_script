#!/usr/bin/env python3

import sys
import argparse
import numpy as np
from ase.io import read, write
from ase.build import make_supercell
from ase.build.tools import sort

# 扩胞
def expand_cell(input_structure, expansion_matrix):
    atoms = read(input_structure)
    
    # 生成扩胞结构
    supercell = make_supercell(atoms, expansion_matrix)
    
    # 对生成的超胞进行按元素排序
    sorted_supercell = sort(supercell)
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

    return parser.parse_args()

if __name__ == "__main__":
    info = '''python3 supercell_expansion.py -i POSCAR -ds 2 2 2
-ds 2 2 2 表示对角方向的 2x2x2 扩胞
python3 supercell_expansion.py -i POSCAR -nds 0 1 1 1 0 1 1 1 0
-nds 0 1 1 1 0 1 1 1 0 生成一个 3x3 的非对角扩胞

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

    # 生成扩胞
    supercell = expand_cell(input_structure, supercell_matrix)

    # 输出扩胞后的结构文件
    output_file = 'POSCAR_supercell'
    write(output_file, supercell)
    print(f"Supercell structure written to {output_file}")