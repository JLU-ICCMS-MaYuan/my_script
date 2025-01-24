#!/usr/bin/env python3
import os
import shutil
import argparse
from pathlib import Path

import spglib
from ase.io import read, write
from ase import Atoms

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
    
    
def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Rename and move VASP files.")
    parser.add_argument("-b", "--begin_id", type=int, required=True, 
        help="Starting index for naming")
    parser.add_argument(
        "-w", "--way_of_naming", choices=["retain_old_name", "formula_symmetry",], 
        default="retain_old_name", 
        help="The naming convention to use (default is 'retain_old_name', you can set another of 'formula_symmetry')."
    )
    parser.add_argument('-p', '--prec',  type=float, nargs='?', default=0.01, 
        help="Symmetry precision (default is 0.01, use default precision values)")
    parser.add_argument('-i', '--input-file-or-directory', default='.', type=str, nargs='+',
        help="The structure file (e.g., POSCAR, CONTCAR) or a directory")
    parser.add_argument('-c', '--detailed-conditions', default=['OUTCAR'], type=str, nargs='+', 
        help="List of subdirectories or patterns for detailed selection")
    parser.add_argument('-od', '--preferred-order', nargs='+', type=str, default=None,
                        help='Custom element preferred_order for sorting, e.g., O C H for ordering oxygen first, followed by carbon and hydrogen')
    parser.add_argument('-t', '--type-of-crystal', type=str, default="prim",
                        help='Custom what type of struct will be output, prim, std, none. none represents the action that the program copy original files to new position')
    return parser.parse_args()

def create_directory(directory):
    """Create a directory if it doesn't exist."""
    if not directory.exists():
        os.mkdir(directory)

def get_vasp_files(input_file_or_directory, detailed_conditions):
    allfiles = []
    for dst_fileordir in input_file_or_directory:
        if os.path.isdir(dst_fileordir):
            for root, dirs, files_in_dir in os.walk(dst_fileordir):
                for filename in files_in_dir:
                    filepath = os.path.join(root, filename)
                    if all(pattern in filepath for pattern in detailed_conditions):
                        allfiles.append(Path(filepath))
        elif os.path.isfile(dst_fileordir):
            allfiles.append(Path(dst_fileordir))
    return allfiles

def generate_new_filename(old_filepath, current_id, way_of_naming, symbols, spacegroup):
    """Generate the new filename based on the naming convention."""
    old_file_name = old_filepath.name
    if way_of_naming == 'retain_old_name':
        new_file_name = f"{current_id}.{old_file_name.split('.')[0]}.vasp"
    elif way_of_naming == 'formula_symmetry':
        new_file_name = f"{current_id}.{str(symbols)}_{spacegroup.replace(' (', '_').replace(')', '_').replace('/', '_')}.vasp"
    return new_file_name

def rename_and_copy_files(vaspfiles, begin_id, way_of_naming, indexed_dirs, prec, preferred_order, type_of_crystal):
    """Rename and copy the files to the new directory."""
    for idx, old_filepath in enumerate(vaspfiles):
        atoms = read(old_filepath)
            # 获取空间群对称性
        lattice = atoms.get_cell()
        positions = atoms.get_scaled_positions()
        numbers = atoms.get_atomic_numbers()
        symbols = atoms.symbols
        cell = (lattice, positions, numbers)
        
        # 使用 spglib 获取空间群
        spacegroup = spglib.get_spacegroup(cell, prec)
        current_id = idx + begin_id
        new_file_name = generate_new_filename(old_filepath, current_id, way_of_naming, symbols, spacegroup)
        new_filepath = indexed_dirs.joinpath(new_file_name)
        print(f"{old_filepath.name} -> {new_file_name}")

        if type_of_crystal == 'std':
            lattice, scaled_positions, numbers = spglib.standardize_cell(cell, symprec=prec)
            std_atoms = sort_by_custom_order(Atoms(cell=lattice, scaled_positions=scaled_positions, numbers=numbers), preferred_order=preferred_order)
            write(new_filepath, std_atoms, direct=True)
        elif type_of_crystal == 'prim':
            lattice, scaled_positions, numbers = spglib.find_primitive(cell, symprec=prec)
            std_atoms = sort_by_custom_order(Atoms(cell=lattice, scaled_positions=scaled_positions, numbers=numbers), preferred_order=preferred_order)
            write(new_filepath, std_atoms, direct=True)
        else:
            shutil.copy(old_filepath, new_filepath)


def main():
    # Parse arguments
    args = parse_arguments()
    
    # Get VASP files
    vaspfiles = get_vasp_files(args.input_file_or_directory, args.detailed_conditions)

    # Create the indexed_dirs directory if it doesn't exist
    indexed_dirs = Path.cwd().joinpath("stdlibs")
    create_directory(indexed_dirs)

    # Rename and copy the files
    rename_and_copy_files(vaspfiles, args.begin_id, args.way_of_naming, indexed_dirs, args.prec, args.preferred_order, args.type_of_crystal)

if __name__ == "__main__":
    main()
