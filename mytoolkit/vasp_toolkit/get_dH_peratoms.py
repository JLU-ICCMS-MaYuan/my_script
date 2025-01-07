#!/usr/bin/env python3
import os
import re
import argparse
import subprocess
from ase.io import read

# 设置命令行参数解析器
def parse_arguments():
    parser = argparse.ArgumentParser(description="Get file paths based on input directory and conditions.")
    parser.add_argument('-i', '--input-file-or-directory', type=str, help="The structure file (e.g., POSCAR, CONTCAR) or a directory")
    parser.add_argument('-c', '--detailed-conditions', type=str, nargs='+', default=[], 
                        help="List of subdirectories or patterns for detailed selection")
    return parser.parse_args()

# 设置需要计算的文件类型，假设所有结构文件都为 .vasp 文件
def calculate_enthalpy(outcar_path):

   try:
        result = subprocess.run(f"grep -s enthalpy {outcar_path} | tail -n 1 | awk '{{print $5}}'", 
                                shell=True, capture_output=True, text=True)
        dH = result.stdout.strip()
        return float(dH)
   except Exception:
        return 6666666666666666666


def get_files(input_file_or_directory, detailed_conditions):
    file_paths = []

    # 检查输入的是文件还是目录
    if os.path.isdir(input_file_or_directory):
        # 如果是目录，遍历所有文件
        for root, dirs, files in os.walk(input_file_or_directory):
            for filename in files:
                filepath = os.path.join(root, filename)
                # 如果文件路径符合所有条件
                if all(pattern in filepath for pattern in detailed_conditions):
                    file_paths.append(filepath)
    elif os.path.isfile(input_file_or_directory):
        # 如果是单个文件
        file_paths.append(input_file_or_directory)
    else:
        print(f"The path {input_file_or_directory} is neither a valid file nor a directory.")

    return file_paths

def get_info(f):
    
    try:
        atom = read(f, format="vasp-out")
    except:
        atom = None
        return None
    
    atom.info['enthalpy'] = calculate_enthalpy(f)
    return atom

def main():
    # 解析命令行参数
    args = parse_arguments()
    input_file_or_directory = args.input_file_or_directory
    detailed_conditions = args.detailed_conditions  # 获取详细的筛选条件

    # 获取符合条件的文件路径
    files = get_files(input_file_or_directory, detailed_conditions)

    # 用来存储文件和计算出来的焓值
    enthalpy_results = []

    # 遍历所有文件并计算焓值
    for f in files:
        atom = get_info(f)
        if atom is not None:
            enthalpy_results.append((f, str(atom.symbols), len(atom), atom.info['enthalpy']/len(atom), atom.get_volume()/len(atom)))
        else:
 
            enthalpy_results.append((f, "None", 6666666666666666666, 6666666666666666666, 6666666666666666666))
    # 按照焓值从小到大排序
    sorted_enthalpy_results = sorted(enthalpy_results, key=lambda x: x[3])

    # 输出结果
    print(f"file    symbol    TotalatomsNumber    Enthalpy(eV/atom)    volume")
    for file, symbol, natom, enthalpy, volume in sorted_enthalpy_results:
        print(f"{file:<}    {symbol:>10}    {natom:>3d}    {enthalpy:>.6f}    {volume:>.6f}")

if __name__ == '__main__':
    main()

