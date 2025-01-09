#!/usr/bin/env python3

import os
import subprocess
import argparse

# 提取POTCAR中的元素顺序
def get_potcar_elements(potcar_file):
    result = subprocess.run(['grep', '-oP', 'VRHFIN =\K\w+', potcar_file], stdout=subprocess.PIPE, text=True)
    return result.stdout.replace('\n', ' ').strip()

# 提取POSCAR中的元素顺序（第六行）
def get_poscar_elements(poscar_file):
    with open(poscar_file, 'r') as f:
        lines = f.readlines()
        return ' '.join(lines[5].split()).strip()

# 比较两个元素顺序
def compare_elements(potcar_elements, poscar_elements):
    return potcar_elements == poscar_elements

# 判断文件是否为空
def is_file_not_empty(file_path):
    return os.path.getsize(file_path) > 0

# 遍历目录并对每个子目录中的POTCAR和POSCAR进行比较
def compare_in_directory(root_dir):
    ele_file = open('elements_order.dat', 'w')
    root_dir = os.path.abspath(root_dir)  # 获取绝对路径
    for root, dirs, files in os.walk(root_dir):
        # 判断是否存在POTCAR和POSCAR文件且文件不为空
        potcar_file = os.path.join(root, 'POTCAR')
        poscar_file = os.path.join(root, 'POSCAR')
        
        if 'POTCAR' in files and 'POSCAR' in files:
            if is_file_not_empty(potcar_file) and is_file_not_empty(poscar_file):
                # 提取POTCAR和POSCAR中的元素顺序
                potcar_elements = get_potcar_elements(potcar_file)
                poscar_elements = get_poscar_elements(poscar_file)

                # 比较元素顺序
                if not compare_elements(potcar_elements, poscar_elements):
                    # print(f"{os.path.abspath(root)}: consistent", file=ele_file)
                    print(f"{os.path.abspath(root)}: not consistent", file=ele_file)
            else:
                print(f"{os.path.abspath(root)}: POTCAR or POSCAR is blank")

    ele_file.close()

    
# 主程序
if __name__ == "__main__":
    # 设置命令行参数解析器
    parser = argparse.ArgumentParser(description="递归比较POTCAR和POSCAR中的元素顺序")
    parser.add_argument('-d', '--directories', type=str, nargs='+', help="要遍历的根目录路径，可以指定多个目录，-d 目录1 目录2 ...")
    
    # 获取命令行参数
    args = parser.parse_args()
    
    # 递归比较每个指定的目录中的文件
    for directory in args.directories:
        compare_in_directory(directory)
