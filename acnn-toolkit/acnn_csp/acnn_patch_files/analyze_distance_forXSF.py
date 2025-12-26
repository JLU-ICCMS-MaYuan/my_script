#!/usr/bin/env python3
import sys
import os
import argparse
from ase.io import read
import numpy as np
from io import StringIO

def extract_xsf_structure(file_path):
    """从非标准 XSF 文件中提取纯晶体结构（CRYSTAL-ATOMS-END 片段）"""
    with open(file_path, 'r', encoding='utf-8') as f:
        lines = f.readlines()
    
    structure_lines = []
    in_structure = False
    for line in lines:
        line_stripped = line.strip()
        if line_stripped.startswith('CRYSTAL'):
            in_structure = True
            structure_lines.append(line)
        elif line_stripped.startswith('END') and in_structure:
            structure_lines.append(line)
            break
        elif in_structure:
            if not any(keyword in line_stripped.upper() for keyword in ['VIRIAL', 'ENERGY', 'STRESS', 'FORCE']):
                structure_lines.append(line)
    
    if not structure_lines:
        print(f"警告：文件 {file_path} 中未找到晶体结构片段", file=sys.stderr)
        return None
    
    return StringIO(''.join(structure_lines))

def get_min_distance(file_path):
    """计算单个结构文件的最小原子距离（兼容非标准 XSF）"""
    structure_io = extract_xsf_structure(file_path)
    if structure_io is None:
        return None
    
    try:
        atoms = read(structure_io, format='xsf')  # 强制按 XSF 格式解析
    except Exception as e:
        print(f"警告：解析文件 {file_path} 失败 - {str(e)}", file=sys.stderr)
        return None
    
    distances = atoms.get_all_distances()
    np.fill_diagonal(distances, 1000)
    return np.min(distances)

def main():
    # 简化参数解析：只保留 -i/--input
    parser = argparse.ArgumentParser(description='计算 XSF 文件中原子的最小距离（仅支持 XSF 格式）')
    parser.add_argument('-i', '--input', required=True, 
                        help='输入 XSF 文件路径或包含 XSF 文件的目录路径')
    args = parser.parse_args()

    input_path = args.input
    file_list = []

    if os.path.isfile(input_path):
        # 处理单个文件（默认按 XSF 解析）
        file_list.append(input_path)
    elif os.path.isdir(input_path):
        # 处理目录：遍历所有文件（无需筛选后缀，ASE 会尝试解析，失败会提示）
        for filename in os.listdir(input_path):
            file_path = os.path.join(input_path, filename)
            if os.path.isfile(file_path):
                file_list.append(file_path)
    else:
        print(f"错误：输入路径 {input_path} 不存在", file=sys.stderr)
        sys.exit(1)

    # 计算并输出结果
    for file_path in file_list:
        min_dist = get_min_distance(file_path)
        if min_dist is not None:
            print(f"{file_path}   {min_dist:.6f}")

if __name__ == "__main__":
    main()