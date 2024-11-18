#!/usr/bin/env python3
import argparse
from ase.io.vasp import read_vasp_xdatcar, write_vasp

# 设置命令行参数解析器
parser = argparse.ArgumentParser(description="Extract a specific frame from XDATCAR and save as POSCAR.")
parser.add_argument("-i", "--input_filename", default="XDATCAR", help="Path to the XDATCAR file (default: 'XDATCAR').")
parser.add_argument('-s', "steps", type=int, help="Step number to extract the frame from XDATCAR.")

# 解析命令行参数
args = parser.parse_args()

# 读取指定的步数帧
frames = read_vasp_xdatcar(args.filename, index=args.steps)

# 写入 POSCAR 文件
write_vasp(f"POSCAR_{args.steps}.vasp", frames)