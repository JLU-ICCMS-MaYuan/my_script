#!/usr/bin/env python3
import argparse
import spglib
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
    parser.add_argument('-i', '--input-file',  type=str, required=True, help='Directory or pattern containing the structure files.')
    parser.add_argument('-od', '--preferred-order', nargs='+', type=str, default=None,
                        help='Custom element preferred_order for sorting, e.g., O C H for ordering oxygen first, followed by carbon and hydrogen')
    parser.add_argument('-p', '--prec', default=0.01, type=float, help="get tolerance of symmetry ")
    return parser.parse_args()


# 主函数
def main():
    # 解析命令行参数
    args = parse_arguments()
    input_file = args.input_file
    preferred_order = args.preferred_order
    prec = args.prec
    # 读取每个文件并添加到atoms_list
    traj_frames = read(input_file, format='traj', index=':')
    for idx, atoms in enumerate(traj_frames):
        # 获取原子总数
        total_atoms = len(atoms)
        # 获取空间群对称性
        lattice = atoms.get_cell()
        positions = atoms.get_scaled_positions()
        numbers = atoms.get_atomic_numbers()
        cell = (lattice, positions, numbers)
        spacegroup = spglib.get_spacegroup(cell, prec)
        
        # 打印原子数和空间群对称性
        print("{:<20} {:<20} {:<20}".format(str(atoms.symbols), total_atoms, spacegroup))
        atoms = sort_by_custom_order(atoms, preferred_order=preferred_order)
        write(f"{idx+1}.{str(atoms.symbols)}.vasp", atoms)

# 运行脚本
if __name__ == '__main__':
    main()
    print("For example: make_traj_file.py -i . -c 100 -o ScH.traj -t CONTCAR   -od Sc H")
