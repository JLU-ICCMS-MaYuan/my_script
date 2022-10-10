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

from math import comb
import re
import os
import random
import logging
from itertools import *
from collections import Counter, defaultdict
from pathlib import Path
import signal
from symbol import atom
import time
from typing import *

import numpy as np
from pyxtal import pyxtal
from pyxtal.tolerance import Tol_matrix
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.ase import AseAtomsAdaptor

from ase.formula import Formula
from ase import Atoms

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
                # print('start alarm signal.')
                r = func(*args, **kwargs)
                # print('close alarm signal.')
                signal.alarm(0)  # 关闭闹钟
                return r
            except RuntimeError as e:
                callback()
        return inner
    return wrap

def after_timeout(): 
    print("Time out!")


class split_wyckoffs:
    '''
    generator_main.py -w ./Ar-Ne-H-spg229-500/ -i ./229.ini method -m mode=specifywps

    input.ini的内容: 
    [splitwps]
    spacegroup_number = 229 
    nameofatoms = ["Ar", "Ne", "H"] 
    wyckoffpositions = {
        '2a': False,
        '6b': False,
        '8c': False,
        '12d': False,
        '12e': True,
        '16f': True,
        '24g': True,
        '24h': True,
        '48i': True,
        '48j': True,
        '48k': True,
        '96l': True,
        }
    nonH_upper_limit = '12e'
    H_lower_limit    = '12d'
    sitesoccupiedrange=[[1,2], 
                        [1,2], 
                        [1,3],] 
    popsize=300 
    maxlimit=150
    distancematrix=[[2.014, 1.908, 1.590],
                    [1.908, 1.802, 1.483],
                    [1.590, 1.483, 1.116],]

    '''
    def __init__(
        self              ,
        work_path: Path,
        spacegroup_number: int,
        nameofatoms: list[str], 

        wyckoffpositions : Dict[str, bool],
        nonH_upper_limit : str,
        H_lower_limit : str,
        sitesoccupiedrange: list[int],
        popsize: int,
        maxlimit: int,
        distancematrix: list[list[float]]
        ):
        
        self.spacegroup_number  = spacegroup_number
        self.nameofatoms        = nameofatoms

        self.wyckoffpositions   = wyckoffpositions
        self.nonH_upper_limit   = nonH_upper_limit
        self.H_lower_limit      = H_lower_limit
        self.sitesoccupiedrange = sitesoccupiedrange

        self.work_path          = Path(work_path)
        if not self.work_path.exists():
            self.work_path.mkdir(parents=True)
        self.popsize            = popsize
        self.maxlimit           = maxlimit
        self.distancematrix     = np.array(distancematrix)
        self._group = self.get_group(
            self.nameofatoms, 
            self.wyckoffpositions,
            self.sitesoccupiedrange,
            self.nonH_upper_limit,
            self.H_lower_limit,
            )
        self.structs = []

    def get_group(
        self,
        nameofatoms,
        wyckoffpositions,
        sitesoccupiedrange,
        nonH_upper_limit,
        H_lower_limit,
        ):

        if   len(nameofatoms) == 2:
            self.binary_hydrides(
                nameofatoms,
                sitesoccupiedrange,
                wyckoffpositions,
                nonH_upper_limit,
                H_lower_limit,
            )
        elif len(nameofatoms) == 3:
            _group = self.ternary_hydrides(
                nameofatoms,
                sitesoccupiedrange,
                wyckoffpositions,
                nonH_upper_limit,
                H_lower_limit,
            )
            return _group
        elif len(nameofatoms) == 4:
            self.quaternary_hydrides()
        

    def binary_hydrides(
        self,
        nameofatoms,
        sitesoccupiedrange,
        wyckoffpositions,
        nonH_upper_limit,
        H_lower_limit,
        ):
        _group = defaultdict(list)
        # 判断 H元素 是否为最后一个元素
        if not nameofatoms[-1] == 'H':
            raise ValueError("you have to set nameofatoms in order as `H element is the last one` !!! Just like nameofatoms=['Ar', 'Ne', 'H']")
        
        nonH1_range = sitesoccupiedrange[:-1]
        H_range = sitesoccupiedrange[-1]
        # 考虑该非氢元素可能占据的wps的个数的范围从 `nonH1_range[0]~nonH1_range[-1]+1`
        # 例如： 
        #   nonH1_range = [1, 3]
        # 那么:
        #   range(nonH1_range[0], nonH1_range[-1]+1) = [1, 2, 3]
        for nonH1_num in range(nonH1_range[0], nonH1_range[-1]+1):
            for nonH1_wps, nonH1_rest in self.get_NonHwps(list(wyckoffpositions.keys()) ,wyckoffpositions, nonH1_num, nonH_upper_limit):

                for H_num in range(H_range[0], H_range[-1]+1):
                    for H_wps, H_rest in self.get_Hwps(list(wyckoffpositions.keys()), nonH1_rest, H_num, H_lower_limit):
                        
                        wyckps  = [nonH1_wps, H_wps]
                        nelems  = self.get_natoms([nonH1_wps, H_wps])
                        formula = ''.join(map(str, chain.from_iterable(zip(nameofatoms, nelems))))
                        formula = Formula(formula).format("metal") # 这里的分子式是没有经过约化的。是多少就是多少，8:8:80 不会变成1:1:10
                        _group[formula].append([nelems, wyckps])
        
        _group = dict(_group)
        return _group

    def ternary_hydrides(
        self,
        nameofatoms,
        sitesoccupiedrange,
        wyckoffpositions,
        nonH_upper_limit,
        H_lower_limit,
        ):
        # 考虑该非氢元素可能占据的wps的个数的范围从 `nonH_list1[0]~nonH_list1[-1]+1`
        # 例如： 
        #   nonH_list1 = [1, 3]
        # 那么:
        #   range(nonH_list1[0], nonH_list1[-1]+1) = [1, 2, 3]
        # 判断 H元素 是否为最后一个元素
        _group = defaultdict(list)

        if not nameofatoms[-1] == 'H':
            raise ValueError("you have to set nameofatoms in order as `H element is the last one` !!! Just like nameofatoms=['Ar', 'Ne', 'H']")
        
        nonH1_range, nonH2_range = sitesoccupiedrange[:-1]
        H_range = sitesoccupiedrange[-1]
        for nonH1_num in range(nonH1_range[0], nonH1_range[-1]+1):
            for nonH1_wps, nonH1_rest in self.get_NonHwps(list(wyckoffpositions.keys()) ,wyckoffpositions, nonH1_num, nonH_upper_limit):

                for nonH2_num in range(nonH2_range[0], nonH2_range[-1]+1):
                    for nonH2_wps, nonH2_rest in self.get_NonHwps(list(wyckoffpositions.keys()), nonH1_rest, nonH2_num, nonH_upper_limit):

                        for H_num in range(H_range[0], H_range[-1]+1):
                            for H_wps, H_rest in self.get_Hwps(list(wyckoffpositions.keys()), nonH2_rest, H_num, H_lower_limit):
                                
                                wyckps  = [nonH1_wps, nonH2_wps, H_wps]
                                nelems  = self.get_natoms([nonH1_wps, nonH2_wps, H_wps])
                                formula = ''.join(map(str, chain.from_iterable(zip(nameofatoms, nelems))))
                                formula = Formula(formula).format("metal") # 这里的分子式是没有经过约化的。是多少就是多少，8:8:80 不会变成1:1:10
                                _group[formula].append([nelems, wyckps])
        
        _group = dict(_group)
        return _group

    def quaternary_hydrides(
        self,
        nameofatoms,
        sitesoccupiedrange,
        wyckoffpositions,
        nonH_upper_limit,
        H_lower_limit,
        ):
        # 考虑该非氢元素可能占据的wps的个数的范围从 `nonH_list1[0]~nonH_list1[-1]+1`
        # 例如： 
        #   nonH_list1 = [1, 3]
        # 那么:
        #   range(nonH_list1[0], nonH_list1[-1]+1) = [1, 2, 3]
        # 判断 H元素 是否为最后一个元素
        _group = defaultdict(list)

        if not nameofatoms[-1] == 'H':
            raise ValueError("you have to set nameofatoms in order as `H element is the last one` !!! Just like nameofatoms=['Ar', 'Ne', 'H']")
        
        nonH1_range, nonH2_range, nonH3_range = sitesoccupiedrange[:-1]
        H_range = sitesoccupiedrange[-1]
        for nonH1_num in range(nonH1_range[0], nonH1_range[-1]+1):
            for nonH1_wps, nonH1_rest in self.get_NonHwps(list(wyckoffpositions.keys()) ,wyckoffpositions, nonH1_num, nonH_upper_limit):

                for nonH2_num in range(nonH2_range[0], nonH2_range[-1]+1):
                    for nonH2_wps, nonH2_rest in self.get_NonHwps(list(wyckoffpositions.keys()), nonH1_rest, nonH2_num, nonH_upper_limit):

                        for nonH3_num in range(nonH2_range[0], nonH2_range[-1]+1):
                            for nonH3_wps, nonH3_rest in self.get_NonHwps(list(wyckoffpositions.keys()), nonH2_rest, nonH2_num, nonH_upper_limit):

                                for H_num in range(H_range[0], H_range[-1]+1):
                                    for H_wps, H_rest in self.get_Hwps(list(wyckoffpositions.keys()), nonH3_rest, H_num, H_lower_limit):
                                        
                                        wyckps  = [nonH1_wps, nonH2_wps, nonH3_wps, H_wps]
                                        nelems  = self.get_natoms([nonH1_wps, nonH2_wps, nonH3_wps, H_wps])
                                        formula = ''.join(map(str, chain.from_iterable(zip(nameofatoms, nelems))))
                                        formula = Formula(formula).format("metal") # 这里的分子式是没有经过约化的。是多少就是多少，8:8:80 不会变成1:1:10
                                        _group[formula].append([nelems, wyckps])
        
        _group = dict(_group)
        return _group

    def get_NonHwps(self, original_wps:List, wps:dict, num:int, nonH_upper_limit:str):
        allwps = list(wps.keys())
        for comb_list in combinations_with_replacement(allwps, num):
            for item in set(comb_list):
                if (wps[item] == False) and (comb_list.count(item) != 1): # False 代表该wps占位只能占据一次
                    break
                if original_wps.index(item) > original_wps.index(nonH_upper_limit):
                    break
            else:
                comb_list = sorted(list(comb_list))
                
                restwps = set(allwps) - set(comb_list)
                rest_dict = {key : value for key, value in wps.items() if key in restwps}
                yield comb_list, rest_dict

    def get_Hwps(self, original_wps:List, wps:dict, num:int, H_lower_limit):
        allwps = list(wps.keys())
        for comb_list in combinations_with_replacement(allwps, num):
            for item in set(comb_list):
                if (wps[item] == False) and (comb_list.count(item) != 1): # False 代表该wps占位只能占据一次
                    break
                if original_wps.index(item) < original_wps.index(H_lower_limit):
                    break
            else:
                comb_list = sorted(list(comb_list))
                
                restwps = set(allwps) - set(comb_list)
                rest_dict = {key : value for key, value in wps.items() if key in restwps}
                
                yield comb_list, rest_dict

    def get_natoms(self, wyckps:List[List[str]]):
        amount_for_every_element = []
        for wyck in wyckps:
            natoms = []
            for wp_multi in wyck:
                natoms.extend(re.findall(r"\d+", wp_multi))
            natoms = list(map(int, natoms))
            amount_for_every_element.append(sum(natoms))
        return amount_for_every_element

    @set_timeout(10, after_timeout)
    def _gen_randomly(self):
        '''
        Function:
            create a structure by randomly choosing a formula from `self._group` dictionary !!!
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
            # logger.info("Now the program will try to create ")
            # logger.info(f"every atoms amouts are {amounts}")
            # logger.info(f"wyckoff positions are {wyck}\n")
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
            logger.info(f"The structure input infomation {amounts} {wyck} can't create a reasonable structure!!!")

    @set_timeout(30, after_timeout)
    def _gen_specify_symbols(self, input_atoms: Atoms):

        spacegroup_number = self.spacegroup_number
        nameofatoms = self.nameofatoms

        fomula = str(input_atoms.symbols.get_chemical_formula('metal'))
        species_amounts_sites = random.choice(self._group[fomula])

        amounts = species_amounts_sites[0]
        wyck = species_amounts_sites[1]

        species_radius = list(zip(nameofatoms, self.distancematrix.diagonal()/2.0))

        tm = Tol_matrix()
        for ele_r1, ele_r2 in combinations(species_radius, 2):
            tm.set_tol(ele_r1[0], ele_r2[0], ele_r1[1]+ele_r2[1])

        struc = pyxtal()
        try:
            # logger.info("Now the program will try to create ")
            # logger.info(f"every atoms amouts are {amounts}")
            # logger.info(f"wyckoff positions are {wyck}\n")
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
            logger.info(f"The {fomula} generated by `_gen_specify_symbols` failed ")
            return None

    @classmethod
    def init_from_config(cls, config_d: dict):

        self = cls(
            work_path=config_d["work_path"],
            spacegroup_number=config_d["spacegroup_number"],
            nameofatoms=config_d["nameofatoms"], 
            wyckoffpositions=config_d["wyckoffpositions"], 
            nonH_upper_limit=config_d["nonh_upper_limit"],
            H_lower_limit=config_d["h_lower_limit"],
            sitesoccupiedrange=config_d["sitesoccupiedrange"],
            popsize=config_d["popsize"],
            maxlimit=config_d["maxlimit"],
            distancematrix=config_d["distancematrix"],
        )
        return self


if __name__ == "__main__":

    spl_wps = split_wyckoffs(
        work_path='./test',
        spacegroup_number=229,
        nameofatoms=['Ar', 'Ne', 'H'], 
        wyckoffpositions = {
        '2a': False,
        '6b': False,
        '8c': False,
        '12d': False,
        '12e': True,
        '16f': True,
        '24g': True,
        '24h': True,
        '48i': True,
        '48j': True,
        '48k': True,
        '96l': True,
        },
        nonH_upper_limit = '12d', # 表示非氢元素上限可以取到 12d
        H_lower_limit    = '12d', # 表示非氢元素下限可以取到 12d
        sitesoccupiedrange=[[1,2], 
                            [1,2],
                            [1,3],],
        popsize=300,
        maxlimit=150,
        distancematrix=[[2.014, 1.908, 1.590],
                        [1.908, 1.802, 1.483],
                        [1.590, 1.483, 1.116],],
    )