#!/public/home/mayuan/miniconda3/envs/cage/bin/python3

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

class CrystalSpecifyWyckoffs:

    def get_H(self, H_occupied_wps, h_lower, h_upper):
        if h_upper < h_lower:
            h_lower, h_upper = h_upper, h_lower
        hydrogen_wps = []
        for i in range(h_lower, h_upper + 1):
            for e in product(H_occupied_wps, repeat=i):
                e = sorted(list(e))
                res = Counter(e)
                if (np.array(list(res.values())) < 2).all():
                    if not hydrogen_wps:
                        hydrogen_wps.append(e)
                    else:
                        for h_wp in hydrogen_wps:
                            h_wp = sorted(list(h_wp))
                            if h_wp == e:
                                break
                        else:
                            hydrogen_wps.append(e)

        hydrogen_atoms_number = []
        for wps_strings in hydrogen_wps:
            atoms = []
            for wp_multi in wps_strings:
                atoms.extend(re.findall(r"\d+", wp_multi))
            atoms = list(map(int, atoms))
            hydrogen_atoms_number.append(sum(atoms))

        return hydrogen_wps, hydrogen_atoms_number

    def get_M(self, X_occupied_wps, X_lower, X_upper):
        if X_upper < X_lower:
            X_lower, X_upper = X_upper, X_lower
        X_wps = []
        for i in range(X_lower, X_upper + 1):
            for e in combinations(X_occupied_wps, i):
                e = sorted(list(e))
                X_wps.append(e)
        X_atoms_number = []
        for wps_strings in X_wps:
            atoms = []
            for wp_multi in wps_strings:
                atoms.extend(re.findall(r"\d+", wp_multi))
            atoms = list(map(int, atoms))
            X_atoms_number.append(sum(atoms))
        return X_wps, X_atoms_number

    def get_X(self, X, X_occupied_wps, X_lower, X_upper):
        if X == 'H':
            return self.get_H(X_occupied_wps, X_lower, X_upper)
        else:
            return self.get_M(X_occupied_wps, X_lower, X_upper)


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


