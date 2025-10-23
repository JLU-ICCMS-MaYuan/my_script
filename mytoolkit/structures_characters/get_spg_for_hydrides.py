#!/usr/bin/env python3
import os
import argparse
import spglib
from ase.io import read, write
from ase import Atoms


def parse_arguments():
    """Set up command-line interface."""
    parser = argparse.ArgumentParser(
        description="Analyze space group of a structure using spglib (optionally removing hydrogen atoms)."
    )
    parser.add_argument(
        "-i",
        "--input-file-or-directory",
        type=str,
        help="The structure file (e.g., POSCAR, CONTCAR) or a directory",
    )
    parser.add_argument(
        "-p",
        "--prec",
        type=float,
        nargs="?",
        default=None,
        help="Symmetry precision (default is None, use default precision values)",
    )
    parser.add_argument(
        "-c",
        "--detailed-conditions",
        type=str,
        nargs="+",
        default=[],
        help="List of subdirectories or patterns for detailed selection",
    )
    parser.add_argument(
        "-od",
        "--preferred-order",
        nargs="+",
        type=str,
        default=None,
        help=(
            "Custom element preferred_order for sorting, e.g., O C H "
            "for ordering oxygen first, followed by carbon and hydrogen"
        ),
    )
    parser.add_argument(
        "--only-metal",
        action="store_true",
        help="Remove H atoms from the in-memory structure before symmetry analysis.",
    )
    return parser.parse_args()


def sort_by_custom_order(atoms, preferred_order=None):
    """Sort atoms according to a custom preferred order of elements."""
    if preferred_order:
        symbols = atoms.get_chemical_symbols()
        old_indices = [[idx, preferred_order.index(symbol)] for idx, symbol in enumerate(symbols)]
        new_indices = sorted(old_indices, key=lambda x: x[1])
        final_indices = [idx for idx, _ in new_indices]
        atomscopy = atoms[final_indices].copy()
        return atomscopy
    tags = atoms.get_chemical_symbols()
    deco = sorted([(tag, i) for i, tag in enumerate(tags)])
    indices = [i for tag, i in deco]
    atomscopy = atoms[indices].copy()
    return atomscopy


def remove_hydrogen_atoms(atoms):
    """Return a copy of atoms without hydrogen atoms."""
    symbols = atoms.get_chemical_symbols()
    keep_indices = [idx for idx, symbol in enumerate(symbols) if symbol != "H"]
    if not keep_indices:
        raise ValueError("All atoms are hydrogen; cannot perform --only-metal analysis.")
    return atoms[keep_indices].copy()


def sanitize_atoms(atoms, only_metal, preferred_order=None):
    """Apply optional hydrogen removal and sorting to the atoms object."""
    if only_metal:
        atoms = remove_hydrogen_atoms(atoms)
    return sort_by_custom_order(atoms, preferred_order=preferred_order)


def get_atom_info(atoms, prec, preferred_order=None):
    """Return primitive and standardized cells along with the space group."""
    lattice = atoms.get_cell()
    positions = atoms.get_scaled_positions()
    numbers = atoms.get_atomic_numbers()
    cell = (lattice, positions, numbers)

    spacegroup = spglib.get_spacegroup(cell, prec)

    lattice, scaled_positions, numbers = spglib.find_primitive(cell, symprec=prec)
    pri_atoms = sort_by_custom_order(
        Atoms(cell=lattice, scaled_positions=scaled_positions, numbers=numbers),
        preferred_order=preferred_order,
    )

    lattice, scaled_positions, numbers = spglib.standardize_cell(cell, symprec=prec)
    std_atoms = sort_by_custom_order(
        Atoms(cell=lattice, scaled_positions=scaled_positions, numbers=numbers),
        preferred_order=preferred_order,
    )
    return pri_atoms, std_atoms, spacegroup


def process_files(input_file_or_directory, detailed_conditions, prec=None, preferred_order=None, only_metal=False):
    """Read files/directories, optionally remove hydrogen, and print space-group info."""
    if os.path.isdir(input_file_or_directory):
        prim_dir = "prim"
        std_dir = "std"
        if not os.path.exists(prim_dir) and prec is not None:
            os.makedirs(prim_dir)
        if not os.path.exists(std_dir) and prec is not None:
            os.makedirs(std_dir)

        num = 1
        for root, _, files_in_dir in os.walk(input_file_or_directory):
            for filename in files_in_dir:
                filepath = os.path.join(root, filename)
                if all(pattern in filepath for pattern in detailed_conditions):
                    try:
                        atoms = read(filepath)
                        atoms = sanitize_atoms(atoms, only_metal=only_metal, preferred_order=preferred_order)
                        if prec:
                            pri_atoms, std_atoms, spacegroup = get_atom_info(atoms, prec, preferred_order=preferred_order)
                            symbols = str(atoms.symbols)
                            print(f"{symbols:<10} {len(atoms):<3}   {spacegroup:<15}")
                            if prec is not None:
                                sg_tag = spacegroup.replace(" (", "_").replace(")", "_").replace("/", "_")
                                new_filename = os.path.join(prim_dir, f"{num}.{symbols}_{sg_tag}.vasp")
                                write(new_filename, pri_atoms, direct=True)
                                new_filename = os.path.join(std_dir, f"{num}.{symbols}_{sg_tag}.vasp")
                                write(new_filename, std_atoms, direct=True)
                                num += 1
                        else:
                            prec_symmetry = []
                            for pval in [1e-1, 1e-2, 1e-3, 1e-5, 1e-9]:
                                pri_atoms, std_atoms, spacegroup = get_atom_info(atoms, pval, preferred_order=preferred_order)
                                prec_symmetry.append(spacegroup)
                            symbols = str(atoms.symbols)
                            print(
                                "{:<10} {:<3}   {:<15}   {:<15}   {:<15}   {:<15}   {:<15}   {}".format(
                                    symbols,
                                    len(atoms),
                                    prec_symmetry[0],
                                    prec_symmetry[1],
                                    prec_symmetry[2],
                                    prec_symmetry[3],
                                    prec_symmetry[4],
                                    filepath,
                                )
                            )
                    except Exception as e:
                        print(f"Error processing file: {filepath} - {e}")

    elif os.path.isfile(input_file_or_directory):
        try:
            atoms = read(input_file_or_directory)
            atoms = sanitize_atoms(atoms, only_metal=only_metal, preferred_order=preferred_order)
            symbols = str(atoms.symbols)
            if prec:
                pri_atoms, std_atoms, spacegroup = get_atom_info(atoms, prec, preferred_order=preferred_order)
                sg_tag = spacegroup.replace(" (", "_").replace(")", "_").replace("/", "_")
                new_filename = os.path.join(f"prim_{symbols}_{sg_tag}.vasp")
                write(new_filename, pri_atoms, direct=True)
                new_filename = os.path.join(f"std_{symbols}_{sg_tag}.vasp")
                write(new_filename, std_atoms, direct=True)
                print(f"{symbols:<10} {len(atoms):<3}   {spacegroup:<10}")
            else:
                prec_symmetry = []
                for pval in [1e-1, 1e-2, 1e-3, 1e-5, 1e-9]:
                    pri_atoms, std_atoms, spacegroup = get_atom_info(atoms, pval, preferred_order=preferred_order)
                    prec_symmetry.append(spacegroup)
                print(
                    "{:<10} {:<3}   {:<15}   {:<15}   {:<15}   {:<15}   {:<15}".format(
                        symbols,
                        len(atoms),
                        prec_symmetry[0],
                        prec_symmetry[1],
                        prec_symmetry[2],
                        prec_symmetry[3],
                        prec_symmetry[4],
                    )
                )
        except Exception as e:
            print(f"Error processing file: {input_file_or_directory} - {e}")

    else:
        print(f"The path {input_file_or_directory} is neither a valid file nor a directory.")


def main():
    args = parse_arguments()
    input_file_or_directory = args.input_file_or_directory
    prec = args.prec
    detailed_conditions = args.detailed_conditions
    preferred_order = args.preferred_order
    only_metal = args.only_metal

    if prec:
        print("{:<10} {:<3}   {:<15}".format("symbols", "Num", str(prec)))
    else:
        print(
            "{:<10} {:<3}   {:<15}   {:<15}   {:<15}   {:<15}   {:<15}".format(
                "symbols", "Num", "1e-1", "1e-2", "1e-3", "1e-5", "1e-9"
            )
        )

    try:
        process_files(
            input_file_or_directory,
            detailed_conditions,
            prec,
            preferred_order,
            only_metal=only_metal,
        )
    except Exception as e:
        print(f"An unexpected error occurred: {e}")


if __name__ == "__main__":
    main()

