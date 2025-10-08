#!/usr/bin/env python

import sys
import numpy as np
import os
import re
import argparse

def read_rcore_from_potcar(element):
    """从 POTCAR-<element> 文件中读取 RCORE 值"""
    potcar_file = f"POTCAR-{element}"
    if not os.path.exists(potcar_file):
        print(f"\nError: RCORE value for '{element}' not provided and '{potcar_file}' not found.")
        sys.exit(1)

    try:
        with open(potcar_file, 'r') as f:
            for line in f:
                # 使用正则表达式查找 RCORE，例如 "RCORE =      1.500"
                match = re.search(r'RCORE\s*=\s*(\d+\.\d+)', line)
                if match:
                    rcore_value = float(match.group(1))
                    print(f"Info: Found RCORE={rcore_value} for element '{element}' in '{potcar_file}'")
                    return rcore_value
    except Exception as e:
        print(f"\nError: Could not read or parse '{potcar_file}': {e}")
        sys.exit(1)
    
    print(f"\nError: 'RCORE' value not found in '{potcar_file}' for element '{element}'.")
    sys.exit(1)

def parse_input(input_string):
    """解析输入字符串，返回元素名称列表和 RCORE 值列表"""
    element_names = []
    rcore_values = []

    print("Info: Parsing input as element-based format.")
    parts = input_string.strip().split()
    if not parts:
        return [], []

    for part in parts:
        if ':' in part:
            # 'Al:1.5' 格式
            try:
                name, value_str = part.split(':', 1)
                element_names.append(name)
                rcore_values.append(float(value_str))
            except ValueError:
                print(f"\nError: Invalid format for '{part}'. Expected 'ElementName:Value'.")
                sys.exit(1)
        else:
            # 'Al' 格式，需要从 POTCAR 读取
            name = part
            element_names.append(name)
            rcore_values.append(read_rcore_from_potcar(name))
            
    return element_names, rcore_values

def main():
    """主函数"""
    parser = argparse.ArgumentParser(
        description="This script calculates a distance matrix based on RCORE values.\nThe final distance unit is Angstrom.",
        epilog="""Input format examples:

1. Element-RCORE pairs (RCORE in Bohr):
   distancematrix.py Al:2.8 Be:1.8 H:1.1

2. Element names only (RCORE read from POTCAR files):
   distancematrix.py Al Be H
   (This requires POTCAR-Al, POTCAR-Be, POTCAR-H files in the current directory)

3. You can also mix formats:
   distancematrix.py Al:2.8 Be H
""",
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument(
        'elements',
        metavar='ELEMENT',
        nargs='+',
        type=str,
        help="One or more element specifications (e.g., Al, Be:1.8, H)"
    )

    args = parser.parse_args()
    # 将接收到的多个参数（一个列表）合并成一个字符串，以兼容 parse_input 函数
    input_arg = " ".join(args.elements)
    print(f"Input arguments combined: \"{input_arg}\"\n")

    # 解析输入参数
    element_names, rcore_values_bohr = parse_input(input_arg)

    if not rcore_values_bohr:
        print("Error: Could not obtain RCORE values from input. Input might be empty.")
        sys.exit(1)

    # --- 计算 ---
    # bohr -> Angstrom, 并乘以0.7的缩放因子
    rcore_angstrom = np.array([rcore_values_bohr]) * 0.529177 * 0.7 
    # 计算距离矩阵
    distance_matrix = rcore_angstrom + rcore_angstrom.T

    # --- 输出 ---
    print("\n-------------------------- OUTPUT --------------------------")
    # 风格1: 保留已有的矩阵输出
    print("Style 1: Distance Matrix (Angstrom)")
    print(distance_matrix)
    print("")

    # 风格2: 新增的键值对输出
    print("Style 2: Pairwise Bond-Length Format (Angstrom)")
    output_pairs = []
    num_elements = len(element_names)
    # 遍历矩阵的上三角部分
    for i in range(num_elements):
        for j in range(i, num_elements):
            elem1 = element_names[i]
            elem2 = element_names[j]
            dist = distance_matrix[i, j]
            output_pairs.append(f"{elem1}-{elem2}={dist:.5f}")
    
    print(",".join(output_pairs))
    print("----------------------------------------------------------")

if __name__ == "__main__":
    main()

