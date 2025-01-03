#!/usr/bin/env python3
import os
import argparse
import spglib
from ase.io import read, write
from ase import Atoms


# 设置命令行参数解析器
def parse_arguments():
    parser = argparse.ArgumentParser(description="Analyze space group of a structure using spglib")
    parser.add_argument('-i', '--input-file-or-directory', type=str, help="The structure file (e.g., POSCAR, CONTCAR) or a directory")
    parser.add_argument('-p', '--prec',  type=float, nargs='?', default=None, 
                        help="Symmetry precision (default is None, use default precision values)")
    parser.add_argument('-c', '--detailed-conditions', type=str, nargs='+', default=[], 
                        help="List of subdirectories or patterns for detailed selection")
    parser.add_argument('-od', '--preferred-order', nargs='+', type=str, default=None,
                        help='Custom element preferred_order for sorting, e.g., O C H for ordering oxygen first, followed by carbon and hydrogen')

    return parser.parse_args()

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

def get_atom_info(atoms, prec, preferred_order=None):
    
    # 获取空间群对称性
    lattice = atoms.get_cell()
    positions = atoms.get_scaled_positions()
    numbers = atoms.get_atomic_numbers()
    cell = (lattice, positions, numbers)
    
    # 使用 spglib 获取空间群
    spacegroup = spglib.get_spacegroup(cell, prec)

    # 获取标准化结构和原始结构
    lattice, scaled_positions, numbers = spglib.find_primitive(cell, symprec=prec)
    pri_atoms = sort_by_custom_order(Atoms(cell=lattice, scaled_positions=scaled_positions, numbers=numbers), preferred_order=preferred_order)
    lattice, scaled_positions, numbers = spglib.standardize_cell(cell, symprec=prec)
    std_atoms = sort_by_custom_order(Atoms(cell=lattice, scaled_positions=scaled_positions, numbers=numbers), preferred_order=preferred_order)
    return pri_atoms, std_atoms, spacegroup

def process_files(input_file_or_directory, detailed_conditions, prec=1e-2, preferred_order=None):

    if os.path.isdir(input_file_or_directory):
        
        prim_dir = "prim"
        std_dir  = "std"
        if not os.path.exists(prim_dir):
            os.makedirs(prim_dir)
        if not os.path.exists(std_dir):
            os.makedirs(std_dir)
        # 如果是目录，遍历所有文件
        
        num = 1
        for root, dirs, files_in_dir in os.walk(input_file_or_directory):
            for filename in files_in_dir:
                filepath = os.path.join(root, filename)
                # print(filepath, pattern)
                # print(all(pattern in filepath for pattern in detailed_conditions))
                if all(pattern in filepath for pattern in detailed_conditions):
                    atoms = sort_by_custom_order(read(filepath), preferred_order=preferred_order)
                    if prec:
                        # 如果指定了精度，使用指定的精度进行分析
                        pri_atoms, std_atoms, spacegroup = get_atom_info(atoms, prec, preferred_order=preferred_order)
                        print("{:<10} {:<3}   {:<15}".format(str(atoms.symbols), len(atoms), spacegroup))
                        new_filename = os.path.join(prim_dir, f"{num}.{str(atoms.symbols)}_{spacegroup.replace(' (', '_').replace(')', '_').replace('/', '_')}.vasp"); write(new_filename, pri_atoms)
                        new_filename = os.path.join(std_dir, f"{num}.{str(atoms.symbols)}_{spacegroup.replace(' (', '_').replace(')', '_').replace('/', '_')}.vasp"); write(new_filename, std_atoms)
                        num += 1
                    else:
                        # 如果没有指定精度，使用默认精度值列表进行分析
                        prec_symmetry = []
                        for p in [1e-1, 1e-2, 1e-3, 1e-5, 1e-9]:
                            pri_atoms, std_atoms, spacegroup = get_atom_info(atoms, p, preferred_order=preferred_order)
                            prec_symmetry.append(spacegroup)
                        print("{:<10} {:<3}   {:<15}   {:<15}   {:<15}   {:<15}   {:<15}  {}".format(str(atoms.symbols), len(atoms), prec_symmetry[0], prec_symmetry[1], prec_symmetry[2], prec_symmetry[3], prec_symmetry[4], filepath))
    
    elif os.path.isfile(input_file_or_directory):
        # 如果是文件，直接处理该文件
        atoms = sort_by_custom_order(read(input_file_or_directory), preferred_order=preferred_order)
        if prec:
            # 如果指定了精度，使用指定的精度进行分析
            pri_atoms, std_atoms, spacegroup = get_atom_info(atoms, prec, preferred_order=preferred_order)
            new_filename = os.path.join(f"prim_{str(atoms.symbols)}_{spacegroup.replace(' (', '_').replace(')', '_').replace('/', '_')}.vasp"); write(new_filename, pri_atoms)
            new_filename = os.path.join(f"std_{str(atoms.symbols)}_{spacegroup.replace(' (', '_').replace(')', '_').replace('/', '_')}.vasp"); write(new_filename, std_atoms)
            print("{:<10} {:<3}   {:<10}".format(str(atoms.symbols), len(atoms), spacegroup))
        else:
            # 如果没有指定精度，使用默认精度值列表进行分析
            prec_symmetry = []
            for p in [1e-1, 1e-2, 1e-3, 1e-5, 1e-9]:
                pri_atoms, std_atoms, spacegroup = get_atom_info(atoms, p, preferred_order=preferred_order)
                prec_symmetry.append(spacegroup)
            print("{:<10} {:<3}   {:<15}   {:<15}   {:<15}   {:<15}   {:<15}".format(str(atoms.symbols), len(atoms), prec_symmetry[0], prec_symmetry[1], prec_symmetry[2], prec_symmetry[3], prec_symmetry[4]))
                
    else:
        print(f"The path {input_file_or_directory} is neither a valid file nor a directory.")

def main():
    # 解析命令行参数
    args = parse_arguments()
    input_file_or_directory = args.input_file_or_directory
    prec = args.prec
    detailed_conditions = args.detailed_conditions  # 获取详细的筛选条件
    preferred_order = args.preferred_order
    
    if prec:
        print("{:<10} {:<3}   {:<15}".format("symbols", "Num", str(prec)))
    else:
        print("{:<10} {:<3}   {:<15}   {:<15}   {:<15}   {:<15}   {:<15}".format("symbols", "Num", "1e-1", "1e-2", "1e-3", "1e-5", "1e-9"))
        
    try:
        process_files(input_file_or_directory, detailed_conditions, prec, preferred_order)
    except Exception as e:
        print(f"Error occurred: {e}")

if __name__ == '__main__':
    main()
