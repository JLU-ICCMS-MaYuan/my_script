#!/usr/bin/env python3
import os
import argparse
from ase.io import read, write
import glob

# 设置命令行参数解析器
def parse_arguments():
    parser = argparse.ArgumentParser(description='Convert multiple structure files to a single traj file.')
    parser.add_argument('-i', '--input_dir', type=str, required=True, help='Directory containing the structure files.')
    parser.add_argument('-o', '--output_file', type=str, required=True, help='Output traj file name.')
    return parser.parse_args()

# 主函数
def main():
    # 解析命令行参数
    args = parse_arguments()
    input_dir = args.input_dir
    output_file = args.output_file

    # 使用glob模块来匹配多种文件格式
    # 匹配 POSCAR, .vasp 和 CONTCAR 文件
    poscar_files = glob.glob(os.path.join(input_dir, 'POSCAR*')) + \
                   glob.glob(os.path.join(input_dir, '*.vasp')) + \
                   glob.glob(os.path.join(input_dir, 'CONTCAR*')) + \
                   glob.glob(os.path.join(input_dir, 'OUTCAR*'))

    # 创建一个列表来存储所有的结构
    atoms_list = []

    # 读取每个文件并添加到atoms_list
    for poscar in poscar_files:
        atoms = read(poscar)
        atoms_list.append(atoms)

    # 将所有结构写入一个traj文件
    write(output_file, atoms_list)
    print(f"All structure files have been saved to {output_file}")

# 运行脚本
if __name__ == '__main__':
    main()
