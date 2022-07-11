#!/work/home/may/miniconda3/bin/python3

"""
使用帮助：
    例子：
    enumeration.py -s 201 -xwp 2a 4b -xpc 1 -hwp 6d 8e 12f 12g 24h -hpc 2 3 -e
        -s 201                  :   space group is 201
        -xwp 2a 2b              :   X element occupied 2a 2b wyckoff positons
        -xpc                    :   at least get 1 wps to permutation and combination
        -hwp 6d 8e 12f 12g 24h  :   H element occupied 6d 8e 12f 12g 24h wyckoff positons
        -hpc 2 3                :   at least get 2 wps to permutation and combination, 
                                    at most get 3 wps to permutation and combination

"""


from pyxtal import pyxtal
# from pyxtal.lattice import Lattice
from pyxtal.symmetry import Group
from pymatgen.io.vasp import Poscar

import re
import os
import itertools as it

from argparse import ArgumentParser


def get_H(
        H_occupied_wps,
        H_permutation_combination_lower_limit,
        H_permutation_combination_upper_limit
        ):
    """
    获得氢原子的所有wp的占位的组合情况 以及 每种组合情况的原子数为多少
    example: 1a, 2b, 3c
        hydrogen_wps          = [['1a'], ['2b'], ['3c'], ['1a', '2b'], ['1a', '3c'], ['2b', '3c']]
        hydrogen_atoms_number = [[ 1  ], [ 2  ], [ 3  ], [ 1,    2  ], [ 1,    3],   [ 2,    3  ]]
    Notes: 不考虑氢占据所有wp这种组合情况，也就是说，不考虑["1a", "1b", "1c"]这种情况
    Returns:  hydrogen_wps
              hydrogen_atoms_number
    """
    if H_permutation_combination_upper_limit is None:
        H_permutation_combination_upper_limit = len(H_occupied_wps)

    hydrogen_wps = []
    for i in range(H_permutation_combination_lower_limit,
                   H_permutation_combination_upper_limit+1):
        for e in it.combinations(H_occupied_wps, i):
            hydrogen_wps.append(list(e))
    hydrogen_atoms_number = []
    for wps_strings in hydrogen_wps:
        atoms = []
        for wp_multi in wps_strings:
            atoms.extend(re.findall("\d+", wp_multi))
        atoms = list(map(int, atoms))
        hydrogen_atoms_number.append(sum(atoms))

    return hydrogen_wps, hydrogen_atoms_number


def get_X(X_occupied_wps, X_permutation_combination):
    """
    获得氢原子的所有wp的占位的组合情况 以及 每种组合情况的原子数为多少
    example: 1a, 2b, 3c
        other_wps      = [['1a'], ['2b'], ['3c'], ['1a', '2b'], ['1a', '3c'], ['2b', '3c']]
        X_atoms_number = [[ 1  ], [ 2  ], [ 3  ], [ 1,    2  ], [ 1,    3],   [ 2,    3  ]]
    Notes: 不考虑氢占据所有wp这种组合情况，也就是说，不考虑["1a", "1b", "1c"]这种情况
    Returns:  other_element_wps
              X_atoms_number
    """
    X_wps = []
    for i in range(X_permutation_combination, len(X_occupied_wps)+1):
        for e in it.combinations(X_occupied_wps, i):
            X_wps.append(list(e))
    X_atoms_number = []
    for wps_strings in X_wps:
        atoms = []
        for wp_multi in wps_strings:
            atoms.extend(re.findall("\d+", wp_multi))
        atoms = list(map(int, atoms))
        X_atoms_number.append(sum(atoms))
    return X_wps, X_atoms_number


# def select_X_element(spg_number, X_occupied_numbers):
#     """
#     input:   输入一个数字 X_occupied_numbers，表示选定 该空间群的多重度按照从低到高顺序的 X_occupied_numbers个wp 被 X元素占据
#     Returns: 返回两个列表, 第一个列表是X元素占据wp的列表，第二个列表是H元素占据wp的列表
#     """
#     gp = Group(spg_number)
#     all_wps_list = gp.get_wp_list(reverse=True)
#     wp_occupied_by_x = all_wps_list[:X_occupied_numbers]
#     wp_occupied_by_h = all_wps_list[X_occupied_numbers:]
#     return wp_occupied_by_x, wp_occupied_by_h


def creat_struct(spacegroup_number, H_wps, h_atoms_number, X_wps, X_atoms_number):

    destination = os.path.join(str(spacegroup_number),  "struc")
    if not os.path.exists(destination):
        os.makedirs(destination)

    struc = pyxtal()
    for wp_x, atoms_x_num in zip(X_wps, X_atoms_number):
        x_species = ["Ca"]
        for wp_h, atoms_h_num in zip(H_wps, h_atoms_number):
            h_species = ["H"]
            all_species = x_species + h_species
            all_atoms = [atoms_x_num, atoms_h_num]
            all_sites = [wp_x] + [wp_h]
            # print(all_species, all_atoms, all_sites)
            # my_lat = Lattice(ltype="Cubic", volume=27, PBC=[1, 1, 1])
            struc.from_random(3,
                              spacegroup_number,
                              all_species,
                              all_atoms,
                              factor=1.0,
                              sites=all_sites,
                              # lattice=my_lat
                              )

            struct_pymatgen = struc.to_pymatgen()

            all_wp_name = ["X"] + wp_x + ["H"] + wp_h
            wp_name = "_".join(all_wp_name)

            struct_name = struct_pymatgen.composition.formula.replace(' ', '') + "---" + wp_name
            print(struct_name)
            filepath = os.path.join(destination, struct_name + ".vasp")
            Poscar(struct_pymatgen).write_file(filepath)


if __name__ == "__main__":

    parser = ArgumentParser()
    parser.add_argument(
        "-s",
        "--SpaceGroupNumber",
        default=None,
        type=int,
        help="请输入一个空间群号, 例如: 45, 93, 230"
    )
    parser.add_argument(
        "-xwp",
        "--x-wyckoff-positions",
        dest="X_occupied_wps",
        action="store",
        default=None,
        type=str,
        nargs='+',
        help="请输入非氢元素占据wps"
    )
    parser.add_argument(
        "-xpc",
        "--xwps-permutation_combination",
        dest="X_permutation_combination",
        default=1,
        type=int,
        help="对所有 非氢元素的wps 进行排列组合"
    )
    parser.add_argument(
        "-hwp",
        "--hydrogen-wyckoff-positions",
        dest="H_occupied_wps",
        action="store",
        default=None,
        type=str,
        nargs='+',
        help="请输入氢元素占据wps"
    )
    parser.add_argument(
        "-hpc",
        "--hwps-permutation-combination",
        dest="H_permutation_combination",
        default=[None, None],
        nargs='+',
        action="store",
        type=int,
        help="对所有 氢元素的wps 进行排列组合, 至少选取hpcl个进行排列组合，至多选取hpcu个进行排列组合"
    )
    parser.add_argument(
        "-e",
        "--enumeration",
        action="store_true",
        default=False,
        help="是否枚举并产生所有结构"
    )
    parser.add_argument(
        "-w",
        "--WriteWorkLog",
        action="store_true",
        default=False,
        help="是否写文件日志"
    )

    args = parser.parse_args()

    spacegroup_number         = args.SpaceGroupNumber
    X_occupied_wps            = args.X_occupied_wps
    X_permutation_combination = args.X_permutation_combination
    H_occupied_wps            = args.H_occupied_wps
    hpc_lower, hpc_upper      = args.H_permutation_combination
    WriteWorkLog              = args.WriteWorkLog
    enumeration               = args.enumeration
    # wps_occupied_by_X, wp_occupied_by_H = select_X_element(spacegroup_number, X_occupied_wp_number)
    X_wps, X_atoms_number = get_X(X_occupied_wps, X_permutation_combination)
    H_wps, h_atoms_number = get_H(H_occupied_wps, hpc_lower, hpc_upper)


    if enumeration:
        creat_struct(spacegroup_number, H_wps, h_atoms_number, X_wps, X_atoms_number)

    if WriteWorkLog:
        if os.path.exists(str(spacegroup_number)):
            path = os.path.abspath(str(spacegroup_number))
            work_log_path = os.path.join(path, "work.log")
            with open(work_log_path, "a") as f:
                f.write("space group number = %s                 \n" % spacegroup_number)
                f.write("X element occupies wps which is \n  %s  \n" % X_wps)
                f.write("The whole X occupation situation is %d  \n" % len(X_wps))
                f.write("H element occupies wps which is \n  %s  \n" % H_wps)
                f.write("The whole H occupation situation is %d  \n" % len(H_wps))
                total = len(X_wps) * len(H_wps)
                f.write("The total structures number is %s     \n\n" % total)


