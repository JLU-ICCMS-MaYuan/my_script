#!/usr/bin/env python3
import argparse
import random

def replace_nb_in_poscar_grouped(poscar_path:str, output_path:str, proportions:dict[str:float]):
    with open(poscar_path, 'r') as file:
        lines = file.readlines()
    print(proportions)
    atom_count = list(map(int, lines[6].strip().split()))
    elements = [element for element in proportions.keys()]
    element_counts = [int(ratio * sum(atom_count)) for ratio in proportions.values()]
    print(atom_count), print(element_counts)
    assert sum(element_counts) == sum(atom_count), "Total atoms number isn't consitant !"

    element_coords = {element: [] for element in elements}
    atom_start_line = 8
    atom_coordinates = lines[atom_start_line:atom_start_line + sum(atom_count)]
    indices = list(range(sum(atom_count)))
    random.shuffle(indices)

    index = 0
    for element, count in zip(elements, element_counts):
        for _ in range(count):
            element_coords[element].append(atom_coordinates[indices[index]])
            index += 1

    lines[5] = ' '.join(elements) + '\n'
    lines[6] = ' '.join(map(str, element_counts)) + '\n'
    new_atom_lines = []
    for element in elements:
        new_atom_lines.extend(element_coords[element])

    lines[atom_start_line:atom_start_line + sum(atom_count)] = new_atom_lines

    with open(output_path, 'w') as file:
        file.writelines(lines)

def main():
    parser = argparse.ArgumentParser(description="Replace elements in POSCAR file according to specified proportions.")
    parser.add_argument("-i", "--inputfile", help="Directory of the POSCAR files.")
    parser.add_argument("-o", "--outputfile", help="Directory to save the output files.")
    parser.add_argument("-e", "--elements", nargs='+', help="List of elements to be used in the alloy.")
    parser.add_argument("-p", "--proportion", nargs='+', help="Proportion of each element in the alloy (e.g., 0.25 for equal proportions).")
    
    args = parser.parse_args()

    # Initialize proportions dictionary
    proportions = {el: float(prop) for el, prop in zip(args.elements, args.proportion)}

    replace_nb_in_poscar_grouped(args.inputfile, args.outputfile, proportions)

if __name__ == "__main__":
    main()
