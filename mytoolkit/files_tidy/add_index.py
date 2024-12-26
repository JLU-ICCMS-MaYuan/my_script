#!/usr/bin/env python3
import os
import shutil
from pathlib import Path
from ase.io import read
import spglib
import argparse

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Rename and move VASP files.")
    parser.add_argument("-b", "--begin_id", type=int, required=True, help="Starting index for naming")
    parser.add_argument(
        "-w", "--way_of_naming", choices=["retain_old_name", "formula_symmetry",], 
        default="retain_old_name", 
        help="The naming convention to use (default is 'retain_old_name', you can set another of 'formula_symmetry')."
    )
    parser.add_argument('-p', '--prec',  type=float, nargs='?', default=0.01, 
        help="Symmetry precision (default is 0.01, use default precision values)")
    return parser.parse_args()

def create_directory(directory):
    """Create a directory if it doesn't exist."""
    if not directory.exists():
        os.mkdir(directory)

def get_vasp_files():
    """Get all VASP files in the current directory."""
    vaspfiles = list(Path.cwd().glob("*.vasp"))
    return sorted(vaspfiles)

def generate_new_filename(old_filepath, current_id, way_of_naming, symbols, spacegroup):
    """Generate the new filename based on the naming convention."""
    old_file_name = old_filepath.name
    if way_of_naming == 'retain_old_name':
        new_file_name = f"{current_id}.{old_file_name.split('.')[0]}.vasp"
    elif way_of_naming == 'formula_symmetry':
        new_file_name = f"{current_id}.{str(symbols)}_{spacegroup.replace(' (', '_').replace(')', '_').replace('/', '_')}.vasp"
    return new_file_name

def rename_and_copy_files(vaspfiles, begin_id, way_of_naming, indexed_dirs, prec):
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
        shutil.copy(old_filepath, new_filepath)

def main():
    # Parse arguments
    args = parse_arguments()
    
    # Get VASP files
    vaspfiles = get_vasp_files()

    # Create the indexed_dirs directory if it doesn't exist
    indexed_dirs = Path.cwd().joinpath("stdlibs")
    create_directory(indexed_dirs)

    # Rename and copy the files
    rename_and_copy_files(vaspfiles, args.begin_id, args.way_of_naming, indexed_dirs, args.prec)

if __name__ == "__main__":
    main()
