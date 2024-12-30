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

    return parser.parse_args()

def get_atom_info(atoms, prec):
    total_atoms = len(atoms)
    # 获取空间群对称性
    lattice = atoms.get_cell()
    positions = atoms.get_scaled_positions()
    numbers = atoms.get_atomic_numbers()
    cell = (lattice, positions, numbers)
    
    # 使用 spglib 获取空间群
    spacegroup = spglib.get_spacegroup(cell, prec)

    # 获取标准化结构和原始结构
    lattice, scaled_positions, numbers = spglib.find_primitive(cell, symprec=prec)
    pri_atoms = Atoms(cell=lattice, scaled_positions=scaled_positions, numbers=numbers)
    
    lattice, scaled_positions, numbers = spglib.standardize_cell(cell, symprec=prec)
    std_atoms = Atoms(cell=lattice, scaled_positions=scaled_positions, numbers=numbers)
    
    return pri_atoms, std_atoms, spacegroup

def process_files(input_file_or_directory, detailed_conditions, prec):

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
                    atoms = read(filepath)
                    if prec:
                        # 如果指定了精度，使用指定的精度进行分析
                        pri_atoms, std_atoms, spacegroup = get_atom_info(atoms, prec)
                        print("{:<10} {:<3}   {:<15}".format(str(atoms.symbols), len(atoms), spacegroup))
                        new_filename = os.path.join(prim_dir, f"{num}.{str(atoms.symbols)}_{spacegroup.replace(' (', '_').replace(')', '_').replace('/', '_')}.vasp"); write(new_filename, pri_atoms)
                        new_filename = os.path.join(std_dir, f"{num}.{str(atoms.symbols)}_{spacegroup.replace(' (', '_').replace(')', '_').replace('/', '_')}.vasp"); write(new_filename, std_atoms)
                        num += 1
                    else:
                        # 如果没有指定精度，使用默认精度值列表进行分析
                        prec_symmetry = []
                        for p in [1e-1, 1e-2, 1e-3, 1e-5, 1e-9]:
                            pri_atoms, std_atoms, spacegroup = get_atom_info(atoms, p)
                            prec_symmetry.append(spacegroup)
                        print("{:<10} {:<3}   {:<15}   {:<15}   {:<15}   {:<15}   {:<15}  {}".format(str(atoms.symbols), len(atoms), prec_symmetry[0], prec_symmetry[1], prec_symmetry[2], prec_symmetry[3], prec_symmetry[4], filepath))
    
    elif os.path.isfile(input_file_or_directory):
        # 如果是文件，直接处理该文件
        atoms = read(input_file_or_directory)
        
        if prec:
            # 如果指定了精度，使用指定的精度进行分析
            pri_atoms, std_atoms, spacegroup = get_atom_info(atoms, prec)
            new_filename = os.path.join(f"prim_{str(atoms.symbols)}_{spacegroup.replace(' (', '_').replace(')', '_').replace('/', '_')}.vasp"); write(new_filename, pri_atoms)
            new_filename = os.path.join(f"std_{str(atoms.symbols)}_{spacegroup.replace(' (', '_').replace(')', '_').replace('/', '_')}.vasp"); write(new_filename, std_atoms)
            print("{:<10} {:<3}   {:<10}".format(str(atoms.symbols), len(atoms), spacegroup))
        else:
            # 如果没有指定精度，使用默认精度值列表进行分析
            prec_symmetry = []
            for p in [1e-1, 1e-2, 1e-3, 1e-5, 1e-9]:
                pri_atoms, std_atoms, spacegroup = get_atom_info(atoms, p)
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

    if prec:
        print("{:<10} {:<3}   {:<15}".format("symbols", "Num", str(prec)))
    else:
        print("{:<10} {:<3}   {:<15}   {:<15}   {:<15}   {:<15}   {:<15}".format("symbols", "Num", "1e-1", "1e-2", "1e-3", "1e-5", "1e-9"))
        
    try:
        process_files(input_file_or_directory, detailed_conditions, prec)
    except Exception as e:
        print(f"Error occurred: {e}")

if __name__ == '__main__':
    main()
