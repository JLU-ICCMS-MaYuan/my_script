#!/usr/bin/env python3
import numpy as np
import argparse
from pymatgen.io.lobster.outputs import Cohpcar
from pymatgen.electronic_structure.core import Orbital, Spin


def average_cohp_for_orbitals(num_bonds, energy_grid, input_orbit_name1, input_orbit_name2):
    total_cohp = np.zeros_like(energy_grid)
    total_icohp = np.zeros_like(energy_grid)
    icohp_atEF  = 0.0
    fermi_index = np.argmin(np.abs(energy_grid - 0))
    fermi_energy = energy_grid[fermi_index]
    
    total_orbit_num = 0
    for atom_pair_idx, orbitals_dict in cohpcar.orb_res_cohp.items(): # orbitals_dict = {1:{}, 2:{}, ...., n:{}}, 1,2,3,...是原子对的编号
        total_orbit_num = 0
        for orb_name, lobster_data in orbitals_dict.items(): # lobster_data = {5s-1s:{}, 6s-1s:{}, ...., 4f_7-1s:{}}
            orb_name1 = orb_name.split('-')[0]
            orb_name2 = orb_name.split('-')[1]
            if input_orbit_name1 in orb_name1 and input_orbit_name2 in orb_name2: 
                total_orbit_num += 1
                # lobster_data = {'COHP':{}, 'ICOHP':{}, 'orbitals':{}, 'length':{}, 'sites':{}}
                for spin_type, cohp in lobster_data['COHP'].items(): 
                # lobster_data['COHP'] = {'Spin.up:1':{...}, 'Spin.down:2':{...}}
                # spin_type如果考虑了自旋，就有Spin.up和Spin.down的区分，如果不考虑自旋，就只有Spin.up
                    total_cohp += cohp  # 加总
                for spin_type, icohp in lobster_data['ICOHP'].items(): 
                    total_icohp += icohp
        else:
            # 遍历完编号为atom_pair_idx的原子对之后，检查找到的指定的轨道数的数量
            if total_orbit_num == 0: # 如果为0， 说明指定的轨道名称可能有问题
                print(f"In atom_pair {atom_pair_idx}, there are {total_orbit_num} {input_orbit_name1}-{input_orbit_name2}. Something wrong")
    else:
        total_cohp = total_cohp/num_bonds
        total_icohp = total_icohp/num_bonds
        icohp_atEF = total_icohp[fermi_index]
        print(f"orbital pair of {input_orbit_name1}-{input_orbit_name2}: {total_orbit_num}, ICOHP at the Fermi level: {icohp_atEF:2.3f}")
    return total_cohp, total_icohp, icohp_atEF



if __name__ == "__main__":
    info = """Sum COHP for given orbital pairs.
Such as: For the system of CeSc2H24


analysiscohp.py -o '5s-1s' '6s-1s' '5py-1s' '5pz-1s' '5px-1s' '5dxy-1s' '5dyz-1s' '5dz2-1s' '5dxz-1s' '5dx2-1s' '4f_3-1s' '4f_2-1s' '4f_1-1s' '4f0-1s' '4f1-1s' '4f2-1s' '4f3-1s'
analysiscohp.py -o 5s-1s 6s-1s 5p-1s 5d-1s 4f-1s


num_bonds: 30
num_orbis: 17
name_orbis: ['5s-1s', '6s-1s', '5py-1s', '5pz-1s', '5px-1s', '5dxy-1s', '5dyz-1s', '5dz2-1s', '5dxz-1s', '5dx2-1s', '4f_3-1s', '4f_2-1s', '4f_1-1s', '4f0-1s', '4f1-1s', '4f2-1s', '4f3-1s']
orbital pair of 5s-1s: 1, ICOHP at the Fermi level: -0.057,        
orbital pair of 6s-1s: 1, ICOHP at the Fermi level: -0.104,        
orbital pair of 5p-1s: 3, ICOHP at the Fermi level: -0.018,        
orbital pair of 5d-1s: 5, ICOHP at the Fermi level: -0.801,        
orbital pair of 4f-1s: 7, ICOHP at the Fermi level: -0.150,        
average ICOHP of the pair at the Fermi level: -1.130,
If you can plus all projected orbital pair to check the results.
    """
    parser = argparse.ArgumentParser(description=info)
    parser.add_argument("-f", "--file", type=str, default="COHPCAR.lobster", help="Path to COHPCAR.lobster file")
    parser.add_argument("-o", "--orbitals", nargs="+", required=True,
                        help="List of orbital pairs to analyze, e.g., 1s-5s 1s-6s 1s-5p 1s-5d 1s-4f")

    args = parser.parse_args()
    
    cohpcar    = Cohpcar(filename="COHPCAR.lobster")
    num_bonds  = len(cohpcar.cohp_data)-1;print(f"num_bonds: {num_bonds}")
    num_orbis  = len(cohpcar.orb_res_cohp['1'].values());print(f"num_orbis: {num_orbis}")
    name_orbis = list(cohpcar.orb_res_cohp['1'].keys());print(f"name_orbis: {name_orbis}")
    energy_grid   = cohpcar.energies
    
    proj_cohp_matrix = [energy_grid]  # 第一列是能量
    proj_cohp_label = ["energy"]  # 对应的列名
    for orb_pair in args.orbitals:
        try:
            input_orbit_name1, input_orbit_name2 = orb_pair.split("-")
        except ValueError:
            print(f"Invalid orbital pair format: {orb_pair}, expected format like 4s-4f or 4f_1-1s")
            continue

        total_cohp, total_icohp, icohp_atEF = average_cohp_for_orbitals(num_bonds, energy_grid, input_orbit_name1, input_orbit_name2)
        proj_cohp_matrix.extend([total_cohp, total_icohp])  # 添加到列中
        proj_cohp_label.extend([orb_pair+'-COHP', orb_pair+'-ICOHP'])
    proj_cohp_matrix = np.array(proj_cohp_matrix).T # 不转置，每一行是一个轨道组的cohp，转置之后，每一列是一个轨道组的cohp
    with open("proj_cohp.dat", "w") as f:
        # 写标题
        f.write("  ".join(proj_cohp_label) + "\n")
        # 写数据
        for row in proj_cohp_matrix:
            f.write("  ".join(f"{val:12.6f}" for val in row) + "\n")
    
    
    average_cohp  = cohpcar.cohp_data['average']['COHP'][Spin.up]
    average_icohp = cohpcar.cohp_data['average']['ICOHP'][Spin.up]
    fermi_index = np.argmin(np.abs(energy_grid - 0))
    average_icohp_atEF = average_icohp[fermi_index]
    cohp_matrix  = np.array([energy_grid, average_cohp, average_icohp]).T 
    print(f"average ICOHP of the pair at the Fermi level: {average_icohp_atEF:2.3f}")
    cohp_label = ["energy", "COHP", "ICOHP"]
    with open("average_cohp.dat", "w") as f:
        # 写标题
        f.write("  ".join(cohp_label) + "\n")
        # 写数据
        for row in cohp_matrix:
            f.write("  ".join(f"{val:12.6f}" for val in row) + "\n")
    print("If you can plus all projected orbital pair to check the results. The average ICOHP is divided by the total number of atomic pairs, rather than by the number of orbitals.")