#!/usr/bin/env python3
import os
import argparse
from ase.io import read, write

# 设置命令行参数解析器
def parse_arguments():
    parser = argparse.ArgumentParser(description="Merge multiple TRAJ files into a single TRAJ file.")
    parser.add_argument('-i', '--input_files', nargs='+', help="List of input TRAJ files to merge")
    parser.add_argument('-o', '--output', required=True, help="Output file name (e.g., merged.traj)")
    return parser.parse_args()

def main():
    # 解析命令行参数
    args = parse_arguments()
    input_files = args.input_files
    output_file = args.output

    atoms_list = []
    for filename in input_files:
        # 读取每个traj文件的所有帧
        traj_frames = read(filename, format='traj', index=':')
        print('filename:{}\nframes found={}\nframes read={}'.format(filename, len(traj_frames), len(traj_frames)))
        atoms_list.extend(traj_frames)

    print('Total frames: {}'.format(len(atoms_list)))

    # 将所有结构写入一个新的traj文件
    write(output_file, atoms_list)
    print(f"All frames have been saved to {output_file}")

if __name__ == '__main__':
    main()
