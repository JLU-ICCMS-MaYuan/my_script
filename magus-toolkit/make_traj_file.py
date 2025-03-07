#!/usr/bin/env python3
import os
import glob
import argparse
from ase.io import read, write

def sort_by_custom_order(atoms, preferred_order=None):
    """
    Sort atoms according to a custom preferred preferred_order of elements.

    Parameters:
    atoms: ASE Atoms object
        The atomic structure to be sorted.
    preferred_order: list of str
        The desired preferred_order of elements, e.g., ['O', 'C', 'H'].

    Returns:
    ASE Atoms object
        The sorted atomic structure.
    """

    if preferred_order:
        symbols = atoms.get_chemical_symbols()
        old_indices = [[idx, preferred_order.index(symbol)] for idx, symbol in enumerate(symbols)]
        new_indices = sorted(old_indices, key=lambda x: x[1])
        final_indices = [idx for idx, symbol_idx in new_indices]
        atomscopy = atoms[final_indices].copy()
        return atomscopy
    else:
        tags = atoms.get_chemical_symbols()
        deco = sorted([(tag, i) for i, tag in enumerate(tags)])
        indices = [i for tag, i in deco]
        atomscopy = atoms[indices].copy()
        return atomscopy

# 设置命令行参数解析器
def parse_arguments():
    parser = argparse.ArgumentParser(
        description=(
            'Convert multiple structure files to a single traj file.\n'
            'Example usage:\n'
            '  process_traj.py -i . -c 200 -o Ba.traj -t POSCAR CONTCAR\n'

        )
    )
    parser.add_argument('-i', '--input-dirs',          type=str, nargs='+', required=True, help='Directory or pattern containing the structure files.')
    parser.add_argument('-c', '--detailed-conditions',type=str, nargs='+', required=True, help='List of subdirectories or patterns for detailed selection')
    parser.add_argument('-o', '--output-file',        type=str, required=True, help='Output traj file name')
    parser.add_argument('-t', '--file-types',          type=str, nargs='+', default=['POSCAR*'], help='List of type of file types, default=[POSCAR*], more situations is *vasp ')
    parser.add_argument('-od', '--preferred-order', nargs='+', type=str, default=None,
                        help='Custom element preferred_order for sorting, e.g., O C H for ordering oxygen first, followed by carbon and hydrogen')
    return parser.parse_args()


def search_files(input_dirs, file_types, detailed_conditions):
    # 搜索所有匹配的文件

    if os.path.exists("dst_files.dat"):
        os.remove("dst_files.dat")
        
    f = open("dst_files.dat", 'a')
    files = []
    for file_type in file_types:
        # print(glob.glob(os.path.join(input_dirs, '**', file_type), recursive=True))
        for input_dir in input_dirs:
            for dsttype_file in glob.glob(os.path.join(input_dir, '**', file_type), recursive=True):
                for condition in detailed_conditions:
                    if condition not in dsttype_file:
                        break
                else:
                    files.append(dsttype_file)
                    f.write(os.path.abspath(dsttype_file)+'\n')
    f.close()

    return files


# 主函数
def main():
    # 解析命令行参数
    args = parse_arguments()
    input_dirs = args.input_dirs
    output_file = args.output_file
    detailed_conditions = args.detailed_conditions
    file_types = args.file_types
    preferred_order = args.preferred_order
    
    # 使用glob模块来匹配文件类型，支持递归查找
    # 遍历多个用户指定的子目录或模式

    dst_files = search_files(input_dirs, file_types, detailed_conditions)
    
    # 创建一个列表来存储所有的结构
    atoms_list = []

    # 读取每个文件并添加到atoms_list
    for file in dst_files:
        print(file)
        atoms = sort_by_custom_order(read(file), preferred_order=preferred_order)
        atoms_list.append(atoms)

    # 将所有结构写入一个traj文件
    write(output_file, atoms_list)
    print(f"All structure files have been saved to {output_file}")

# 运行脚本
if __name__ == '__main__':
    main()
    print("For example: make_traj_file.py -i . -c 100 -o ScH.traj -t CONTCAR   -od Sc H")
