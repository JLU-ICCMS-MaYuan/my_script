#!/usr/bin/env python3
import numpy as np
import argparse
from pymatgen.io.lobster.outputs import Cohpcar
from pymatgen.electronic_structure.core import Orbital, Spin


def average_cohp_for_orbitals(num_bonds, energy_grid, input_orbit_name1, input_orbit_name2):
    """
    Calculates the average COHP and ICOHP for a given pair of orbitals.
    """
    total_cohp = np.zeros_like(energy_grid)
    total_icohp = np.zeros_like(energy_grid)
    fermi_index = np.argmin(np.abs(energy_grid - 0))
    
    total_orbit_num = 0
    # Iterate over each bond (atom pair) in the COHPCAR
    for atom_pair_idx, orbitals_dict in cohpcar.orb_res_cohp.items():
        # Iterate over each orbital pair contribution for the current bond
        for orb_name, lobster_data in orbitals_dict.items():
            orb_name1, orb_name2 = orb_name.split('-')
            # Check if the orbital names match the user's input
            if input_orbit_name1 in orb_name1 and input_orbit_name2 in orb_name2: 
                total_orbit_num += 1
                # Sum up COHP and ICOHP values
                for spin_type, cohp in lobster_data['COHP'].items(): 
                    total_cohp += cohp
                for spin_type, icohp in lobster_data['ICOHP'].items(): 
                    total_icohp += icohp
    
    if total_orbit_num > 0:
        # Average the values over the number of bonds
        total_cohp /= num_bonds
        total_icohp /= num_bonds
        icohp_atEF = total_icohp[fermi_index]
    else:
        print(f"Warning: No contributions found for orbital pair {input_orbit_name1}-{input_orbit_name2}")
        icohp_atEF = 0.0

    return total_cohp, total_icohp, icohp_atEF


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run COHP analysis on a COHPCAR.lobster file.", 
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-f", "--file", type=str, default="COHPCAR.lobster", 
                        help="Path to the COHPCAR.lobster file.")
    parser.add_argument("-o", "--orbitals", nargs="+", required=True,
                        help="List of orbital pairs to analyze, such as:\n"
                             "Method 1:\n"
                             "5s-1s 6s-1s 5p-1s 5d-1s 4f-1s\n"
                             "5s-3s 5s-3p 5s-3d 5s-4s 6s-3s 6s-3p 6s-3d 6s-4s 5p-3s 5p-3p 5p-3d 5p-4s 5d-3s 5d-3p 5d-3d 5d-4s4f-3s 4f-3p 4f-3d 4f-4s\n"
                             "Method 2:\n"
                             "5s-1s 6s-1s 5py-1s 5pz-1s 5px-1s 5dxy-1s 5dyz-1s 5dz2-1s 5dxz-1s 5dx2-1s 4f_3-1s 4f_2-1s 4f_1-1s 4f0-1s 4f1-1s 4f2-1s 4f3-1s\n"
                             "3s-1s 3px-1s 3py-1s 3pz-1s 3dxy-1s 3dyz-1s 3dz2-1s 3dxz-1s 3dx2-1s\n"
                        )
    parser.add_argument("-vd", "--Vd", type=float, help="Lower energy bound (Vd) for ICOHP range analysis.")
    parser.add_argument("-vm", "--Vm", type=float, help="Upper energy bound (Vm) for ICOHP range analysis.")

    args = parser.parse_args()
    
    # --- Data Loading and Initial Processing ---
    try:
        cohpcar = Cohpcar(filename=args.file)
    except Exception as e:
        print(f"Error reading {args.file}: {e}")
        exit()
        
    num_bonds = len(cohpcar.cohp_data) - 1
    energy_grid = cohpcar.energies
    print(f"Analysis based on {num_bonds} bonds found in {args.file}")
    
    proj_cohp_matrix = [energy_grid]
    proj_cohp_label = ["energy"]
    icohp_at_EF_list = []
    
    # --- Orbital-Resolved COHP Calculation ---
    for orb_pair in args.orbitals:
        try:
            input_orbit_name1, input_orbit_name2 = orb_pair.split("-")
        except ValueError:
            print(f"Invalid orbital pair format: '{orb_pair}'. Please use a format like '4s-4f' or '4f_1-1s'.")
            continue

        total_cohp, total_icohp, icohp_atEF = average_cohp_for_orbitals(num_bonds, energy_grid, input_orbit_name1, input_orbit_name2)
        proj_cohp_matrix.extend([total_cohp, total_icohp])
        proj_cohp_label.extend([orb_pair + '-COHP', orb_pair + '-ICOHP'])
        icohp_at_EF_list.append(icohp_atEF)

    proj_cohp_matrix = np.array(proj_cohp_matrix).T
    
    with open("proj_cohp.dat", "w") as f:
        f.write("  ".join(proj_cohp_label) + "\n")
        np.savetxt(f, proj_cohp_matrix, fmt='%12.6f')
            
    # --- Average COHP Calculation ---
    average_cohp = cohpcar.cohp_data['average']['COHP'][Spin.up]
    average_icohp = cohpcar.cohp_data['average']['ICOHP'][Spin.up]
    fermi_index = np.argmin(np.abs(energy_grid - 0))
    average_icohp_atEF = average_icohp[fermi_index]
    
    cohp_matrix = np.array([energy_grid, average_cohp, average_icohp]).T
    
    with open("average_cohp.dat", "w") as f:
        f.write("  ".join(["energy", "COHP", "ICOHP"]) + "\n")
        np.savetxt(f, cohp_matrix, fmt='%12.6f')

    # --- Output ICOHP at Fermi Level (EF) ---
    print("\n--- ICOHP values at the Fermi level (E_Fermi = 0.000 eV) ---")
    for i, orb_pair in enumerate(args.orbitals):
        print(f"  {orb_pair+'-ICOHP':<20}: {icohp_at_EF_list[i]:12.6f}")
    print(f"  {'average-ICOHP':<20}: {average_icohp_atEF:12.6f}")
    print("-" * 55)

    # --- Output ICOHP in Specified Energy Range [Vd, Vm] ---
    if args.Vd is not None and args.Vm is not None:
        vd_index = np.argmin(np.abs(energy_grid - args.Vd))
        vm_index = np.argmin(np.abs(energy_grid - args.Vm))
        actual_vd = energy_grid[vd_index]
        actual_vm = energy_grid[vm_index]
        
        print(f"\n--- ICOHP values in energy range [{actual_vd:.3f} eV, {actual_vm:.3f} eV] ---")
        print(f"{'Orbital Pair':<20} | {'ICOHP @ EF':>12} | {'ICOHP @ Vd':>12} | {'ICOHP @ Vm':>12} | {'ICOHP(Vm-Vd)':>14}")
        print("-" * 78)

        # Projected orbital pairs
        for i, orb_pair in enumerate(args.orbitals):
            label = orb_pair + '-ICOHP'
            col_index = proj_cohp_label.index(label)
            
            icohp_at_vd = proj_cohp_matrix[vd_index, col_index]
            icohp_at_vm = proj_cohp_matrix[vm_index, col_index]
            icohp_diff = icohp_at_vm - icohp_at_vd
            
            print(f"{label:<20} | {icohp_at_EF_list[i]:12.6f} | {icohp_at_vd:12.6f} | {icohp_at_vm:12.6f} | {icohp_diff:14.6f}")

        # Average ICOHP
        avg_icohp_at_vd = average_icohp[vd_index]
        avg_icohp_at_vm = average_icohp[vm_index]
        avg_icohp_diff = avg_icohp_at_vm - avg_icohp_at_vd
        print(f"{'average-ICOHP':<20} | {average_icohp_atEF:12.6f} | {avg_icohp_at_vd:12.6f} | {avg_icohp_at_vm:12.6f} | {avg_icohp_diff:14.6f}")
        print("-" * 78)