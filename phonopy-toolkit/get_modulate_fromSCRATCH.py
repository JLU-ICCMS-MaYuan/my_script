#!/usr/bin/env python3
import argparse
import pprint
import os
import shutil
from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms
from phonopy.interface.vasp import read_vasp, write_vasp
from phonopy.file_IO import parse_FORCE_SETS

def parse_matrix_arg(matrix_input: list[float]) -> list[list[float]]:
    """Parse a matrix list into a 3x3 list of floats."""
    if len(matrix_input) != 9:
        raise ValueError(f"Invalid matrix format: Expected 9 numbers, but got {len(matrix_input)}.")

    # Reshape the list into a 3x3 matrix
    matrix = [
        [matrix_input[0], matrix_input[1], matrix_input[2]],
        [matrix_input[3], matrix_input[4], matrix_input[5]],
        [matrix_input[6], matrix_input[7], matrix_input[8]],
    ]
    pprint.pprint(matrix)
    return matrix


def generate_force_constant(
    _unitcell: PhonopyAtoms,
    _supercell_matrix: list[list[float]] = [[2, 0, 0], [0, 2, 0], [0, 0, 2]],
    _primitive_matrix: list[list[float]] = [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
) -> Phonopy:
    phonon = Phonopy(
        unitcell=_unitcell,
        supercell_matrix=_supercell_matrix,
        primitive_matrix=_primitive_matrix,
    )
    force_sets = parse_FORCE_SETS(filename="FORCE_SETS")
    phonon.dataset = force_sets  # 另一种比较早的写法是, phonon.set_displacement_dataset(force_sets), 不过也快被抛弃了
    phonon.produce_force_constants()
    return phonon


def generate_modulation(
    _ph_input: Phonopy,
    _supercell_matrix: list[list[float]] = [[2, 0, 0], [0, 2, 0], [0, 0, 2]],
    _q_point: list[float] = [0.5, 0.5, 0],
    _band_index: int = 2,
    _amplitudes: list[float] = [-0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3],
    _phase: float = 0,
):
    """Calculate modulation.
    
    Parameters
    ----------
    _ph_input : Phonopy
        Phonopy object containing the force constants and other data.
    _supercell_matrix : list[list[float]], optional
        Supercell matrix, by default [[2, 0, 0], [0, 2, 0], [0, 0, 2]].
    _q_point : list[float], optional
        q-point in fractional coordinates, by default [0.5, 0.5, 0].
    _band_index : int, optional
        Band index of the phonon mode, by default 2.
    _amplitudes : list[float], optional
        List of amplitudes for the modulation, by default [0.1, 0.2, 0.3].
    _phase : float, optional
        Phase factor for the modulation, by default 0.
    """
    _phonon_modes = [[_q_point, _band_index, amp, _phase] for amp in _amplitudes]
    
    _ph_input.run_modulations(
        dimension=_supercell_matrix, 
        phonon_modes=_phonon_modes,
    )
    cells = _ph_input.get_modulated_supercells()
    return cells


def write_modulated_supercells(
    modulations: list[PhonopyAtoms], 
    amplitudes: list[float],
    q_point: list[float],
    band_index: int,):
    """Write the modulated supercells to POSCAR files in a 'modulations' directory."""
    # Create the 'modulations' directory if it doesn't exist
    modulations_dir = "modulations"+"_"+"qpoints_"+str(q_point[0])+"_"+str(q_point[1])+"_"+str(q_point[2])+"_"+"band_index_"+str(band_index)
    os.makedirs(modulations_dir, exist_ok=True)
    
    # Write each modulated supercell to a POSCAR file
    i = 1
    for phonoatoms, amp in zip(modulations, amplitudes):
        file_name = f"{i}.modu_delta{amp:.2f}.vasp"
        file_path = os.path.join(modulations_dir, file_name)
        write_vasp(file_path, phonoatoms)
        i += 1


def generate_amplitudes(start: float, stop: float, step: float) -> list[float]:
    """Generate a list of amplitudes based on start, stop, and step."""
    amplitudes = []
    current = start
    while current <= stop:
        amplitudes.append(current)
        current += step
    return amplitudes


def delete_old_modulated_files(
    q_point: list[float],
    band_index: int,
    ):
    """Delete the 'modulations' directory and all its contents."""
    modulations_dir = "modulations"+"_"+"qpoints_"+str(q_point[0])+"_"+str(q_point[1])+"_"+str(q_point[2])+"_"+"band_index_"+str(band_index)
    if os.path.exists(modulations_dir):
        try:
            shutil.rmtree(modulations_dir)
            print(f"Deleted directory and its contents: {modulations_dir}")
        except Exception as e:
            print(f"Error deleting {modulations_dir}: {e}")
    else:
        print(f"Directory does not exist: {modulations_dir}")

if __name__ == "__main__":
    # Set up argparse
    parser = argparse.ArgumentParser(description="Generate force constants and modulation for a given POSCAR.")
    parser.add_argument(
        "-pi", 
        "--poscar-init",
        type=str, 
        help="Path to the POSCAR file."
    )
    parser.add_argument(
        "-smf",
        "--supercell_matrix-for-force_constants",
        nargs="+",
        type=float,
        default=[2.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 2.0],
        help="Supercell matrix as 9 numbers (default: 2 0 0 0 2 0 0 0 2).",
    )
    parser.add_argument(
        "-smm",
        "--supercell_matrix-for-modulation",
        nargs="+",
        type=float,
        default=[2.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 1.0],
        help="Supercell matrix as 9 numbers (default: 2 0 0 0 2 0 0 0 2).",
    )

    parser.add_argument(
        "-pm",
        "--primitive_matrix",
        type=float,
        nargs="+",
        default=[1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0],
        help="Primitive matrix as 9 numbers (default: 1 0 0 0 1 0 0 0 1).",
    )
    parser.add_argument(
        "-qp",
        "--q_point",
        nargs="+",
        type=float,
        default=[0.5, 0.5, 0.0],
        help="q-point in fractional coordinates (default: 0.5 0.5 0).",
    )
    parser.add_argument(
        "-bi",
        "--band_index",
        type=int,
        default=2,
        help="Band index of the phonon mode (default: 2).",
    )
    parser.add_argument(
        "-ph",
        "--phase",
        type=float,
        default=0,
        help="Phase factor for the modulation (default: 0).",
    )
    parser.add_argument(
        "-amp",
        "--amplitudes",
        nargs=3,
        type=float,
        default=[-3.0, 3.0, 0.2],
        help="Three values: start, stop, and step for amplitudes (default: -3.0 3.0 0.2).",
    )
    args = parser.parse_args()

    # Delete old modulated files
    delete_old_modulated_files(q_point=args.q_point, band_index=args.band_index)

    # Generate amplitudes based on start, stop, and step
    amplitudes = generate_amplitudes(args.amplitudes[0], args.amplitudes[1], args.amplitudes[2])

    # Read POSCAR
    unitcell = read_vasp(args.poscar_init)

    # Parse matrices
    print("supercell_matrix_for_force_constants:")
    supercell_matrix_for_force_constants = parse_matrix_arg(args.supercell_matrix_for_force_constants)
    print("primitive_matrix:")
    primitive_matrix = parse_matrix_arg(args.primitive_matrix)
    print("supercell_matrix_for_modulation:")
    supercell_matrix_for_modulation = parse_matrix_arg(args.supercell_matrix_for_modulation)
    print("qpoint:\n{}".format(args.q_point))
    # Generate force constants and modulation
    phonon = generate_force_constant(unitcell, supercell_matrix_for_force_constants, primitive_matrix)
    modulations = generate_modulation(
        phonon, 
        supercell_matrix_for_modulation, 
        args.q_point, 
        args.band_index, 
        amplitudes, 
        args.phase
    )

    # Write the modulated supercells to POSCAR files
    write_modulated_supercells(modulations=modulations, amplitudes=amplitudes, q_point=args.q_point, band_index=args.band_index)
