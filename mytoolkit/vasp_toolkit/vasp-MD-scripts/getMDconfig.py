#!/usr/bin/env python3
import argparse
from ase.io.vasp import read_vasp_xdatcar, write_vasp, write_vasp_xdatcar

print("begin number is 0, Not include end number. For example: \n  getMDconfig.py -i XDATCAR -b 0 -e 9\n Then you can get 9 structures, not 10 structures.")

# 设置命令行参数解析器
parser = argparse.ArgumentParser(description="Extract a specific frame from XDATCAR and save as POSCAR.")
parser.add_argument("-i", "--filename", default="XDATCAR", help="Path to the XDATCAR file (default: 'XDATCAR').")
parser.add_argument("-b", "--begin", type=int, help="Step begin number to extract the frame from XDATCAR.")
parser.add_argument("-e", "--end", type=int, help="Step end number to extract the frame from XDATCAR.")
# 元素排序顺序
parser.add_argument('-od', '--order', nargs='+', type=str, default=None,
                    help='Custom element order for sorting, e.g., O C H for ordering oxygen first, followed by carbon and hydrogen, You have to make sure the order you specified is consistant with the order of input files')

# 解析命令行参数
args = parser.parse_args()

# 读取指定的步数帧
if args.begin == args.end:
    frames = read_vasp_xdatcar(args.filename, index=args.end)
    # 写入 POSCAR 文件
    if args.order:
        # 按照指定的元素顺序排序
        element_to_index = {element: idx for idx, element in enumerate(args.order)}
        # 获取排序后的索引
        sorted_indices = sorted(range(len(frames)), key=lambda i: element_to_index[frames.get_chemical_symbols()[i]])
        # 根据指定顺序排列原子
        frames = frames[sorted_indices]
        write_vasp(f"poscar_{args.begin}_{args.end}", frames)
    else:
        write_vasp(f"poscar_{args.begin}_{args.end}", frames, sort=True)
else:
    frames = read_vasp_xdatcar(args.filename, slice(args.begin, args.end))
    # 写入 POSCAR 文件
    write_vasp_xdatcar(f"xdatcar_{args.begin}_{args.end}", frames)
