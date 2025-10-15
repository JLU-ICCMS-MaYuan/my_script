#!/usr/bin/env python

import argparse
import os
import sys
import numpy as np
import itertools

# --- Data Section for USPEX-like method ---
# This dictionary contains empirical atomic volumes at 0 GPa (in Angstrom^3/atom) for the first 86 elements.
EMPIRICAL_ATOMIC_VOLUMES = {
    # Period 1
    'H': 5.8, 'He': 20.0,
    # Period 2
    'Li': 21.6, 'Be': 8.1, 'B': 7.5, 'C': 5.7, 'N': 9.2, 'O': 6.9, 'F': 5.5, 'Ne': 22.1,
    # Period 3
    'Na': 39.4, 'Mg': 23.2, 'Al': 16.6, 'Si': 20.0, 'P': 19.5, 'S': 23.9, 'Cl': 20.5, 'Ar': 37.0,
    # Period 4
    'K': 75.3, 'Ca': 43.2, 'Sc': 25.0, 'Ti': 17.7, 'V': 13.8, 'Cr': 12.0, 'Mn': 12.1,
    'Fe': 11.7, 'Co': 11.1, 'Ni': 10.9, 'Cu': 11.8, 'Zn': 15.2, 'Ga': 19.5, 'Ge': 22.6,
    'As': 21.4, 'Se': 27.2, 'Br': 29.0, 'Kr': 48.2,
    # Period 5
    'Rb': 92.2, 'Sr': 56.4, 'Y': 33.0, 'Zr': 23.2, 'Nb': 17.9, 'Mo': 15.6, 'Tc': 14.1,
    'Ru': 13.5, 'Rh': 13.7, 'Pd': 14.7, 'Ag': 17.0, 'Cd': 21.5, 'In': 26.1, 'Sn': 27.0,
    'Sb': 30.2, 'Te': 33.8, 'I': 36.9, 'Xe': 59.8,
    # Period 6
    'Cs': 116.3, 'Ba': 63.9,
    # Lanthanides
    'La': 37.3, 'Ce': 34.5, 'Pr': 34.4, 'Nd': 34.2, 'Pm': 33.9, 'Sm': 32.9, 'Eu': 48.0,
    'Gd': 32.4, 'Tb': 31.6, 'Dy': 31.2, 'Ho': 30.8, 'Er': 30.5, 'Tm': 30.1, 'Yb': 41.3, 'Lu': 29.5,
    # Rest of Period 6
    'Hf': 22.2, 'Ta': 18.0, 'W': 15.8, 'Re': 14.7, 'Os': 14.0, 'Ir': 14.1, 'Pt': 15.1,
    'Au': 16.9, 'Hg': 23.0, 'Tl': 28.5, 'Pb': 30.3, 'Bi': 35.3, 'Po': 35.8, 'At': 39.0, 'Rn': 50.5
}

def setup_arg_parser():
    """Sets up the argument parser for the script with short aliases."""
    parser = argparse.ArgumentParser(
        description="Generate CALYPSO input.dat files for a range of compositions.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    # --- Core Arguments ---
    parser.add_argument(
        '-e', '--elements', nargs='+', required=True, type=str,
        help="List of element names.\nExample: -e Si O"
    )
    parser.add_argument(
        '-r', '--radii', nargs='+', required=True, type=float,
        help="List of covalent radii for each element (in Bohr).\nRequired for distance matrix.\nExample: -r 1.11 0.60"
    )
    parser.add_argument(
        '-c', '--composition', nargs='+', required=True, type=str,
        help="Range of atoms for each element (e.g., '1-5', '3').\nExample: -c 1-2 4"
    )
    parser.add_argument(
        '-p', '--popsize', required=True, type=int,
        help="Population size for the CALYPSO run (e.g., 30)."
    )
    # --- Volume Calculation Method Arguments ---
    parser.add_argument(
        '-m', '--method', type=str, choices=['radius', 'uspex'], default='radius',
        help="Method for volume calculation: 'radius' (default) or 'uspex'."
    )
    parser.add_argument(
        '-P', '--pressure', type=float,
        help="Target pressure in GPa (required for '--method uspex').\n(Note: Uppercase 'P' to avoid conflict)"
    )
    return parser

def calculate_volume_radius(radii_bohr, composition):
    """Estimates the cell volume based on covalent radii (original method)."""
    scaled_radii_angstrom = radii_bohr * 0.529177 * 0.7
    sphere_volumes = (4/3) * np.pi * np.power(scaled_radii_angstrom, 3)
    total_volume = np.sum(sphere_volumes * composition)
    return total_volume * 1.1

def calculate_volume_uspex(elements, composition, pressure_gpa):
    """Estimates the cell volume using the USPEX-like empirical method."""
    v0 = 0.0
    for el, num in zip(elements, composition):
        if el not in EMPIRICAL_ATOMIC_VOLUMES:
            raise ValueError(f"Error: Empirical volume for element '{el}' is not in the database.")
        v0 += EMPIRICAL_ATOMIC_VOLUMES[el] * num
    
    if pressure_gpa is None or pressure_gpa == 0:
        return v0
    
    a = 0.08  # Empirical constant for pressure response
    b = 2.5   # Empirical constant for curvature
    volume_p = v0 * (1.0 / (1.0 + a * pressure_gpa))**(1.0/b)
    return volume_p

def calculate_distance_matrix(radii_bohr):
    """Calculates the NxN matrix of interatomic distances from radii."""
    distance_matrix_angstrom = (radii_bohr[:, np.newaxis] + radii_bohr) * 0.529177 * 0.7
    return distance_matrix_angstrom

def generate_input_file(elements, radii, composition, popsize, num_elements, method, pressure):
    """Generates a single input.dat file for a given composition."""
    # --- Volume Calculation ---
    if method == 'radius':
        volume = calculate_volume_radius(radii, composition)
    elif method == 'uspex':
        volume = calculate_volume_uspex(elements, composition, pressure)
    
    # --- Distance Calculation ---
    distance_matrix = calculate_distance_matrix(radii)

    # --- Directory and File Setup ---
    dir_name = "".join([f"{el}{num}" for el, num in zip(elements, composition)])
    os.makedirs(dir_name, exist_ok=True)
    output_path = os.path.join(dir_name, 'input.dat')
    
    # --- File Content Generation ---
    system_name = "-".join(elements)
    atom_names = " ".join(elements)
    atom_counts = " ".join(map(str, composition))
    distance_lines = [" ".join([f"{dist:.2f}" for dist in row]) for row in distance_matrix]

    # --- Writing to input.dat ---
    with open(output_path, 'w') as f:
        f.write(f"SystemName = {system_name}\n")
        f.write(f"NumberOfSpecies = {num_elements}\n")
        f.write(f"NameOfAtoms = {atom_names}\n")
        f.write(f"NumberOfAtoms = {atom_counts}\n")
        f.write("NumberOfFormula = 1 1\n")
        f.write(f"Volume = {volume:.4f}\n\n")
        f.write("@DistanceOfIon\n")
        f.write("\n".join(distance_lines) + "\n")
        f.write("@End\n\n")
        f.write(f"PopSize = {popsize}\n")
        f.write("SpeSpaceGroup = 8 230\n")
        f.write("Split = T\n")
        f.write("VSC = F\n")
        f.write("MaxNumAtom = 80\n")
        f.write("Command = sh submit.sh\n")
        f.write("MaxStep = 10\n")

def main():
    """Main function to parse arguments and generate files."""
    parser = setup_arg_parser()
    args = parser.parse_args()

    # --- Argument Validation ---
    if args.method == 'uspex' and args.pressure is None:
        parser.error("--pressure (-P) is required when using --method uspex")
    
    num_elements = len(args.elements)
    if not (len(args.radii) == num_elements and len(args.composition) == num_elements):
        print("Error: The number of arguments for --elements (-e), --radii (-r), and --composition (-c) must be the same.")
        sys.exit(1)

    # --- Parse Composition Ranges ---
    composition_ranges = []
    try:
        for r_str in args.composition:
            if '-' in r_str:
                start, end = map(int, r_str.split('-'))
                composition_ranges.append(range(start, end + 1))
            else:
                num = int(r_str)
                composition_ranges.append(range(num, num + 1))
    except ValueError as e:
        print(f"Error: Invalid format in composition range '{r_str}': {e}")
        sys.exit(1)

    # --- Generate All Composition Combinations ---
    all_compositions = list(itertools.product(*composition_ranges))
    total_files = len(all_compositions)
    print(f"Found {total_files} composition(s). Generating files using '{args.method}' method...")

    # *** NEW: Initialize a list to store all directory names ***
    all_dir_names = []

    # --- Loop and Generate Files ---
    for i, comp in enumerate(all_compositions):
        current_composition = np.array(comp)
        dir_name = "".join([f"{el}{num}" for el, num in zip(args.elements, current_composition)])
        
        # *** NEW: Add the generated directory name to our list ***
        all_dir_names.append(dir_name)

        print(f"[{i+1}/{total_files}] Generating input for {dir_name}...")
        generate_input_file(
            elements=args.elements,
            radii=np.array(args.radii),
            composition=current_composition,
            popsize=args.popsize,
            num_elements=num_elements,
            method=args.method,
            pressure=args.pressure
        )
    
    # *** NEW: Write all collected directory names to composition.dat ***
    try:
        output_filename = "composition.dat"
        with open(output_filename, 'w') as f:
            for name in all_dir_names:
                f.write(f"./{name}\n")
        print(f"\nSuccessfully created all {total_files} input files.")
        print(f"All {len(all_dir_names)} directory paths have been written to '{output_filename}'.")
    except IOError as e:
        print(f"\nWarning: Could not write to file '{output_filename}': {e}")
        print("Input files were still created successfully.")


if __name__ == "__main__":
    main()
