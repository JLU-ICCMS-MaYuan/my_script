#!/usr/bin/env python3
import argparse
from pathlib import Path
from ase.io import read
from ase.io.wien2k import write_struct
from ase.spacegroup import crystal
from ase.formula import Formula


def convert(input_file_path, dst_format):
    if input_file_path is None or dst_format is None:
        raise ValueError("你必须指定输入文件路径和输出文件格式")

    input_file_path = Path(input_file_path)
    work_path = Path('.')  # 默认输出到当前目录
    atoms = read(input_file_path)
    atoms_primitive = crystal(symbols=atoms, primitive_cell=True)
    fm = atoms.get_chemical_formula("metal")  # 拿到一个formula的字符串，但它可能不是最简单的形式
    reduced_formula = Formula(fm).reduce()[0].format()
    # 对Ca2H12取它的最简单的形式reduce()[0]，返回结果是一个Formula类，再对这个formula类按照metal的方式输出为字符串
    
    if "cif" in dst_format:
        if "cif" == dst_format:
            dst_file = work_path.joinpath(reduced_formula + ".cif")
            atoms_primitive.write(dst_file, format="cif")
        else:
            atoms_primitive.write(dst_format, format="cif")
    elif "vasp" in dst_format:
        if "vasp" == dst_format:
            dst_file = work_path.joinpath(reduced_formula + ".vasp")
            atoms_primitive.write(dst_file, format="vasp")
        else:
            atoms_primitive.write(dst_format, format="vasp")
    elif "struct" in dst_format:
        if "vasp" == dst_format:
            dst_file = work_path.joinpath(reduced_formula + ".struct")
            write_struct(dst_file, atoms_primitive)
        else:
            write_struct(dst_format, atoms_primitive)


def main():
    parser = argparse.ArgumentParser(description="Convert input structure to specified format.")
    
    # 输入文件路径
    parser.add_argument('-i', '--input_file', type=str, help="Path to the input file (e.g., .cif, .vasp, etc.)")
    
    # 输出文件格式
    parser.add_argument('-o', '--output_format', type=str,
                        help="Output file format: 'cif', 'vasp', or 'struct'")
    
    args = parser.parse_args()

    convert(args.input_file, args.output_format)


if __name__ == "__main__":
    main()
