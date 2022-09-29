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
from pathlib import Path
import signal
import time

import numpy as np
from pyxtal import pyxtal
from pyxtal.tolerance import Tol_matrix
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.ase import AseAtomsAdaptor

from ase.formula import Formula

from psolib.utils.sort_atoms import sort_atoms

logger = logging.getLogger(__name__)

def handle(timeout, frame):  # 收到信号 SIGALRM 后的回调函数，第一个参数是信号的数字，第二个参数是the interrupted stack frame.
    print("bad structure")
    raise RuntimeError("run error")

def set_timeout(timeout, callback):
    def wrap(func):
        def inner(*args, **kwargs):
            try:
                signal.signal(signal.SIGALRM, handle)  # 设置信号和回调函数
                signal.alarm(timeout)  # 设置 num 秒的闹钟
                print('start alarm signal.')
                r = func(*args, **kwargs)
                print('close alarm signal.')
                signal.alarm(0)  # 关闭闹钟
                return r
            except RuntimeError as e:
                callback()
        return inner
    return wrap

def after_timeout():  # 超时后的处理函数
    print("Time out!")


class specify_wyckoffs:

    def __init__(
        self              ,
        work_path: Path,
        spacegroup_number: int,
        nameofatoms: list[str], 
        optionalsites: list[list[str]], 
        sitesoccupiedrange: list[int],
        popsize: int,
        maxlimit: int,
        distancematrix: list[list[float]]
        ):
        
        self.spacegroup_number   = spacegroup_number
        self.nameofatoms        = nameofatoms
        self.optionalsites      = optionalsites
        self.sitesoccupiedrange = sitesoccupiedrange
        self.work_path          = Path(work_path)
        if not self.work_path.exists():
            self.work_path.mkdir(parents=True)
        self.popsize            = popsize
        self.maxlimit           = maxlimit
        self.distancematrix     = np.array(distancematrix)
        self._group             = self.get_group(self.optionalsites, self.sitesoccupiedrange)
        self.structs = []

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
    
    @set_timeout(10, after_timeout)
    def _specify_gen(self):
        '''
        create struct by choosing wyckoff positions
        '''
        spacegroup_number = self.spacegroup_number
        nameofatoms = self.nameofatoms
        fomula = random.choice(list(self._group.keys()))
        species_amounts_sites = random.choice(self._group[fomula])

        amounts = species_amounts_sites[0]
        wyck = species_amounts_sites[1]

        species_radius = list(zip(nameofatoms, self.distancematrix.diagonal()/2.0))

        tm = Tol_matrix()
        for ele_r1, ele_r2 in combinations(species_radius, 2):
            tm.set_tol(ele_r1[0], ele_r2[0], ele_r1[1]+ele_r2[1])


        struc = pyxtal()
        try:
            logger.info("Now the program will try to create ")
            logger.info(f"every atoms amouts are {amounts}")
            logger.info(f"wyckoff positions are {wyck}\n")
            struc.from_random(
                3,
                spacegroup_number,
                nameofatoms,
                amounts,
                factor=2.0,
                sites=wyck,
                tm=tm
            )
            struct_pymatgen = struc.to_pymatgen()
            if struct_pymatgen.composition.num_atoms < float(self.maxlimit):
                return struct_pymatgen
        except Exception as e:
            logger.info(f"except {e}")
            logger.info(f"The structure input infomation {amounts} {wyck} can't create a reasonable structure!!!")

    @classmethod
    def init_from_config(cls, config_d: dict):

        self = cls(
            work_path=config_d["work_path"],
            spacegroup_number=config_d["spacegroup_number"],
            nameofatoms=config_d["nameofatoms"], 
            optionalsites=config_d["optionalsites"], 
            sitesoccupiedrange=config_d["sitesoccupiedrange"],
            popsize=config_d["popsize"],
            maxlimit=config_d["maxlimit"],
            distancematrix=config_d["distancematrix"],
        )
        return self
