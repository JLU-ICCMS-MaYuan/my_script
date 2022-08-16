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


import re
import os
import random
import logging
from itertools import product, combinations, chain
from collections import Counter, defaultdict
from argparse import ArgumentParser
from pprint import pprint

import numpy as np
from pyxtal.symmetry import Group
from pyxtal import pyxtal
from pymatgen.io.vasp import Poscar
from ase.formula import Formula

logger = logging.getLogger("GENERATOR")

class CrystalSpecifyWyckoffs:

    def __init__(
        self              ,
        spacegroup_number ,
        nameofatoms       , 
        optionalsites     , 
        sitesoccupiedrange,        
        ):
        
        self.SpaceGroupNumber   = spacegroup_number
        self.nameofatoms        = nameofatoms
        self.optionalsites      = optionalsites
        self.sitesoccupiedrange = sitesoccupiedrange
        self._group             = self.get_group(self.optionalsites, self.sitesoccupiedrange)
        pprint(self._group); input()
        
        for _ in range(50):
            self.__rand_gen()

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

    def get_group(self, optionalsites, sitesoccupiedrange) -> dict:
        """
        Function

        Args:
            optional_sites (List[List[str]]):
                optional wyckoff positions of each elements,
                e.g. [['4a', "6b", "8c"],['12g', '12h', "16g"]],
                     ['12m', '16g'].
            sitesoccupiedrange (List[List[int]]):
                repeat range of each wyckoff position,
                e.g. [[1, 3], [2, 6], [2, 6]].

            Returns:
            dict: {'formula0': [amount, wyck], ...}
        """
        elems = self.nameofatoms  # ['H', 'Ca', 'C']
        _group = defaultdict(list)
        wyck_list, nelems_list = [], []
        for elem, opt_sites, sites_range in zip(
            elems, optionalsites, sitesoccupiedrange
        ):
            wyck, nelems = self.get_X(elem, opt_sites, *sites_range)
            wyck_list.append(
                wyck
            )  # [A, B, C]: [[['4b'], ['4d', '4f']], [['4b'], ['4d', '4f'], ...]
            nelems_list.append(
                nelems
            )  # [A, B, C]: [[4, 8], [4, 8], ...]
        nelems_comb = product(*nelems_list)  # generator
        wyck_comb = product(*wyck_list)
        for nelems, wyck in zip(nelems_comb, wyck_comb):
            formula = ''.join(map(str, chain.from_iterable(zip(elems, nelems))))
            formula = Formula(formula).format("metal")
            _group[formula].append([nelems, wyck])
        _group = dict(_group)
        return _group
    
        def __rand_gen(self):
        '''
        create struct by choosing wyckoff positions
        '''
        spespacegroup = self.SpaceGroupNumber
            name_of_atoms = self.nameofatoms
        fomula = random.choice(list(self._group.keys()))
        species_amounts_sites = random.choice(self._group[fomula])

            amounts = species_amounts_sites[0]
        wyck = species_amounts_sites[1]

        struc = pyxtal()
        for i in range(10):
            logger.info("begin creating")
            try:
                print(
                    name_of_atoms,
                    amounts,
                    wyck,
                )
                struc.from_random(
                    3,
                    spespacegroup,
                    name_of_atoms,
                    amounts,
                    factor=1.0,
                    sites=wyck,
                    # lattice=my_lat
                )
                print("try successfully!")
                _random_atoms = struc.to_ase()
                break
            except Exception as e:
                logger.debug(f"except {e}")
                continue 
        else:
            logger.warning("this structure may not be physical")

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
        "-n",
        "--nameofatoms",
        default=None,
        type=str,
        nargs="+",
        help="请输入元素种类"
    )
    parser.add_argument(
        "-o",
        "--optionalsites",
        default=None,
        type=str,
        nargs="+",
        help="请输入每种元素可以选择的位置"
    )
    parser.add_argument(
        "-r",
        "--sitesoccupiedrange",
        default=None,
        type=str,
        nargs="+",
        help="请输入每种元素占据的wyckoff postion的范围"
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
    spacegroup_number  = args.SpaceGroupNumber
    nameofatoms        = args.nameofatoms
    optionalsites      = args.optionalsites
    sitesoccupiedrange = args.sitesoccupiedrange
    WriteWorkLog       = args.WriteWorkLog
    enumeration        = args.enumeration

    spacegroup_number = 229
    nameofatoms = ["Ca", "Sr", "H"]
    optionalsites = [['2a', '6b'],
                     ['8c','12d'],
                     ['12e', '16f', '24g', '24h', '48i', '48j', '48k']]
    sitesoccupiedrange =[
        [1,2],
        [1,2],
        [1,3],
    ]
    if enumeration: 
        wyck_crys = CrystalSpecifyWyckoffs(
            spacegroup_number ,
            nameofatoms       , 
            optionalsites     , 
            sitesoccupiedrange, 
            )

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


