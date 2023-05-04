#!/usr/bin/env python3
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
import random
import logging
from itertools import product, combinations_with_replacement, chain, combinations
from collections import defaultdict
from pathlib import Path
import signal
from typing import *

import numpy as np
from pyxtal import pyxtal
from pyxtal.tolerance import Tol_matrix
from pyxtal.symmetry import Group

from ase.formula import Formula
from ase import Atoms

from structuregenerator.checkstructure import checkclathrate
from structuregenerator.msg import Formula_Specified_Error

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

def compare(x, y):
    """
    Compare whether two arrays are identical after dividing by their greatest common divisor
    """
    _x = np.array(x)
    _y = np.array(y)

    # 求_x的最大公约数
    xf = 0
    for i, xi in enumerate(_x):
        xf = np.gcd(xf, xi)
    if xf != 0:
        _x = np.array(x) / xf
        
    # 求_y的最大公约数
    yf = 0
    for i, yi in enumerate(_y):
        yf = np.gcd(yf, yi)
    if yf != 0:
        _y = np.array(y) / yf

    if np.allclose(_x, _y):
        return True
    else:
        return False

def append_composition(group, nameofatoms, wyckoffpositons, spacegroup_number, nelem1, nelem2, hydrogen_content_upper_limit):
    '''
    Input Paramaters:
        group: a dictory, its key is formula, its value is a list, which is 
            group[formula] = [spacegroup_number, nelem2, wyckoffpositons]
    
        nameofatoms: the name of element
            nameofatoms = [La, Be, H]

        wyckoffpositons:
            wyckoffpositons = [[2a], [2b, 2c], [6d]]

        spacegroup_number: int
            spacegroup_number = 9
    
        nelem1: User-defined chemical composition typically includes three scenarios
            nelem1 = [1,2,3]
            nelem1 = [1,None,3]
            nelem1 = [None, None, None]

        nelem2: The program arranges Wyckoff positions based on a certain space group's permutation and combination,
                and then output a number of atoms for each element. Usually, it's not the simplest format.
            nelem2 = [2,4,6]
            nelem2 = [7,9,21]  Although the second number is not proportional to the other numbers, 
                               the first number and the third number satisfy a ratio of 1:3
            nelem2 = [2,4,9]

        hydrogen_content_upper_limit: the upper limit of hydrogen content
            hydrogen_content_upper_limit = 0.7
        
    '''
    append_flag = False
    if all(nelem1) == True:
        if compare(nelem1, nelem2):
            formula = ''.join(map(str, chain.from_iterable(zip(nameofatoms, nelem2))))
            formula = Formula(formula).format("metal") # 这里的分子式是没有经过约化的。是多少就是多少，8:8:80 不会变成1:1:10
            current_hydrogen_content = nelem2[-1] / sum(nelem2)
            if current_hydrogen_content >= hydrogen_content_upper_limit:
                group[formula].append([spacegroup_number, nelem2, wyckoffpositons])
                append_flag = True
    elif any(nelem1) == True:
        ids_of_num = [id for id, num in enumerate(nelem1) if num is not None ]
        numberofato_exceptNone = [num for id, num in enumerate(nelem1) if id in ids_of_num]
        nelems_exceptNone = [num for id, num in enumerate(nelem2) if id in ids_of_num]
        if compare(numberofato_exceptNone, nelems_exceptNone):
            formula = ''.join(map(str, chain.from_iterable(zip(nameofatoms, nelem2))))
            formula = Formula(formula).format("metal") # 这里的分子式是没有经过约化的。是多少就是多少，8:8:80 不会变成1:1:10
            current_hydrogen_content = nelem2[-1] / sum(nelem2)
            if current_hydrogen_content >= hydrogen_content_upper_limit:
                group[formula].append([spacegroup_number, nelem2, wyckoffpositons])
                append_flag = True
    else:
        formula = ''.join(map(str, chain.from_iterable(zip(nameofatoms, nelem2))))
        formula = Formula(formula).format("metal") # 这里的分子式是没有经过约化的。是多少就是多少，8:8:80 不会变成1:1:10
        current_hydrogen_content = nelem2[-1] / sum(nelem2)
        if current_hydrogen_content >= hydrogen_content_upper_limit:
            group[formula].append([spacegroup_number, nelem2, wyckoffpositons])
            append_flag = True

    if append_flag:
        return 1
    else:
        return 0

class split_wyckoffs:
    '''
    '''
    def __init__(
        self              ,
        work_path: Path,
        spacegroup_number: List[List[int]],
        nameofatoms: List[str], 
        numberofatoms: List,
        mults_ranges: List[List[int]],
        occupied_number_ranges: List[List[int]],
        popsize: int,
        maxlimit: int,
        distancematrix: List[List[float]],
        clathrate_ratio : float,

        hydrogen_content : float,

        remain_H_ratio_upperstd,
        fraction_of_hydrogen_volume_lowerstd,
        shr_num_avg_lowerstd,
        cage_regularity_avg_upperstd,
        h2h_network_regularity_avg_upperstd,
        h2h_1nearest_lowerstd,
        nonh2nonh_1nearest_lowerstd,
        ):
        
        self.spacegroup_number = spacegroup_number
        self.nameofatoms       = nameofatoms
        self.numberofatoms     = numberofatoms
        self.mults_ranges = mults_ranges
        self.occupied_number_ranges = occupied_number_ranges

        self.work_path          = Path(work_path)
        if not self.work_path.exists():
            self.work_path.mkdir(parents=True)
        self.popsize            = popsize
        self.maxlimit           = maxlimit
        self.distancematrix     = np.array(distancematrix)
        self.clathrate_ratio    = clathrate_ratio
        self.hydrogen_content   = hydrogen_content

        self.remain_H_ratio_upperstd = 0.44
        self.fraction_of_hydrogen_volume_lowerstd = 0.3
        self.shr_num_avg_lowerstd = 1.6
        self.cage_regularity_avg_upperstd = 0.05
        self.h2h_network_regularity_avg_upperstd = 0.15
        self.h2h_1nearest_lowerstd = 0.9
        self.nonh2nonh_1nearest_lowerstd = 2.2

        self.structs = []
        self.satisfied_spgs, self.wyckoffpositions = self.get_spgs(self.spacegroup_number, self.mults_ranges)
        self._group = self.get_group(
            self.satisfied_spgs,
            self.wyckoffpositions,
            self.nameofatoms,
            self.mults_ranges,  
            self.occupied_number_ranges,
            )

    def get_spgs(self, spacegroup_number, mults_ranges):
        '''
        判断给定的一系列空间群spacegroup_number的wyckoff position的多重度是否满足输入指定的非氢元素可占据的上限mults_ranges
        satisfied_spgs:  将满足要求的空间群返回
        wyckoffpositions:  将满足要求的空间群的wyckoff position的占据情况做成字典返回
        '''
        allspgs = [spgnum for spgnums in spacegroup_number for spgnum in range(spgnums[0], spgnums[1]+1)]
        wyckoffpositions = defaultdict(dict)
        unsatified_spgs = []
        for spgnum in allspgs:
            spg = Group(spgnum, dim=3)
            wps_list = spg.get_wp_list()
            for mults in mults_ranges:
                # 拿到最低多重度的wps: wps_list[-1]
                if mults[-1] >= int(re.search(r'\d+', wps_list[-1]).group()):
                    for wpstr in wps_list:
                        wp = spg.get_wyckoff_position(wpstr)
                        # 检查给定的元素的最高多重度 mults[-1] 是否比该空间群下最低多重度还要低
                        if wp.get_dof() > 0:
                            wyckoffpositions[spgnum][wpstr] = True
                        else:
                            wyckoffpositions[spgnum][wpstr] = False
                else:
                    unsatified_spgs.append(spgnum)
                    break

        print("The multiplicity of wyckoff position of these spacegroup numbers is not satisfy with the requirement that `mults_ranges[:,-1]` < the lowest multiplicity of these spacegroups")           
        print(unsatified_spgs)

        satisfied_spgs = [spg for spg in allspgs if spg not in unsatified_spgs]
        print("The spacegroup meeted requirement is")
        print(satisfied_spgs)

        return satisfied_spgs, wyckoffpositions

    def get_group(
        self,
        satisfied_spgs,
        wyckoffpositions,
        nameofatoms,
        mults_ranges,
        occupied_number_ranges,
        ):
            # check whether the file 'composition.json' exist or not
            import json
            composition_json = self.work_path.joinpath("composition.json")
            if composition_json.exists():
                print("\nNote:-----------------------")
                print("    The composition.json file exists, please confirm that you will read it !!!!")
                with open(composition_json, 'r') as f:
                    _group = json.load(f)
                    return _group
            else:
                _group = self.mult_hydirdes(
                    nameofatoms,
                    occupied_number_ranges,
                    mults_ranges,
                    satisfied_spgs,
                    wyckoffpositions,
                )
                with open(composition_json, 'w') as f:
                    json.dump(_group, f)
                return _group

    def mult_hydirdes(
        self,
        nameofatoms,
        occupied_number_range,
        mults_ranges,
        satisfied_spgs,
        wyckoffpositions,

    ):
        _group = defaultdict(list)
        # 判断 H元素 是否为最后一个元素
        if not nameofatoms[-1] == 'H':
            raise ValueError("you have to set nameofatoms in order as `H element is the last one` !!! Just like nameofatoms=['Ar', 'Ne', 'H']")

        all_posibility = 0
        for spg in satisfied_spgs:
            wps_bool = wyckoffpositions[spg]
            wps_name = list(wps_bool.keys())
            total_ele_wps  = []
            for id, ele in enumerate(self.nameofatoms):
                occu_n = occupied_number_range[id]
                mult_n = mults_ranges[id]
                one_ele_wps = self.get_wps(wps_name, wps_bool, occu_n, mult_n)
                total_ele_wps.append(one_ele_wps)

            special_spg_posibility = 0
            for wps in product(*total_ele_wps):
                if self.check_duplicate(wps, wps_bool):
                    nelems  = self.get_natoms(wps)
                    n = append_composition(_group, self.nameofatoms, wps, spg, self.numberofatoms, nelems, self.hydrogen_content)
                    special_spg_posibility = n + special_spg_posibility
            all_posibility = all_posibility + special_spg_posibility
            print(f"when considering the spg-{spg}, {special_spg_posibility} scenarios just are given!")
        print(f"when considering the all spg , {all_posibility} scenarios just are given!")

        _group = dict(_group)
        if _group:
            return _group
        else:  # 有可能指定的配比在当前空间群和wyckoff组合下并不存在，所以需要提示错误
            msg = f"The specified chemical formula {self.nameofatoms} doesn't exist"
            Formula_Specified_Error(msg)

    def check_duplicate(self, wps:List[List[str]], wps_bool:dict):
        flag_wps = list(chain(*wps))
        for wp in flag_wps:
            if wps_bool[wp] == False and (flag_wps.count(wp) > 1):
                return False
        else:
            return True

    def get_wps(self, wps_name:List, wps_bool:dict, occu_num:List[int], mult_n:List[int]):

        floor_mult = mult_n[0]
        upper_mult = mult_n[1]
        combinated_lists = []
        for num in range(occu_num[0], occu_num[1]+1):
            for comb_list in combinations_with_replacement(wps_name, num):
                # print("comb_list", comb_list); input()
                for item in set(comb_list):
                    if (wps_bool[item] == False) and (comb_list.count(item) > 1): # False 代表该wps占位只能占据一次
                        # 检查：在将wps排列组合的列表里面, 检查占据情况为False并且存在1次以上
                        # 那么这个wps排列组合的列表就是不正确的
                        break
                    if (int(re.search(r'\d+', item).group()) < floor_mult) or (int(re.search(r'\d+', item).group()) > upper_mult):
                        # 检查：在将wps排列组合的列表里面, 包含过高 或者or 过低多重度的wps
                        # 那么这个wps排列组合的列表就是不符合要求的
                        break
                else:
                    comb_list = sorted(list(comb_list))
                    combinated_lists.append(comb_list)
        return combinated_lists

    def get_natoms(self, wyckps:List[List[str]]):
        """
        Function:
            get the occu_number of atoms for every element
        """
        amount_for_every_element = []
        for wyck in wyckps:
            natoms = []
            for wp_multi in wyck:
                natoms.extend(re.findall(r"\d+", wp_multi))
            natoms = list(map(int, natoms))
            amount_for_every_element.append(sum(natoms))
        return amount_for_every_element

    # @set_timeout(10, after_timeout)
    def _gen_randomly(self):
        '''
        Function:
            create a structure by randomly choosing a formula from `self._group` dictionary !!!
        ''' 
        formula = random.choice(list(self._group.keys()))
        spg_amounts_sites = random.choice(self._group[formula])
        spg     = spg_amounts_sites[0]
        amounts = spg_amounts_sites[1]
        wyck    = spg_amounts_sites[2]
        nameofatoms = self.nameofatoms
        
        if sum(amounts) > float(self.maxlimit):
            print("Total number is over the limit of max number of atoms")
            print(f"Generating failed !!!")
            return None, None

        # species_radius = list(zip(nameofatoms, self.distancematrix.diagonal()/2.0))
        # tm = Tol_matrix()
        # for ele_r1, ele_r2 in combinations(species_radius, 2):
        #     tm.set_tol(ele_r1[0], ele_r2[0], ele_r1[1]+ele_r2[1])

        struc = pyxtal()
        try:
            print(f"try {spg}-spacegroup {amounts} {wyck}")
            struc.from_random(
                3,
                spg,
                nameofatoms,
                amounts,
                factor=1.5,
                sites=wyck,
            )
            struct_pmg = struc.to_pymatgen()

            if self.clathrate_ratio > np.random.uniform():
                if checkclathrate(
                    struct_pmg,
                    remain_H_ratio_UPPERSTD=self.remain_H_ratio_upperstd,
                    fraction_of_hydrogen_volume_LOWERSTD=self.fraction_of_hydrogen_volume_lowerstd,
                    shr_num_avg_LOWERSTD=self.shr_num_avg_lowerstd,
                    cage_regularity_avg_UPPERSTD=self.cage_regularity_avg_upperstd,
                    h2h_network_regularity_avg_UPPERSTD=self.h2h_network_regularity_avg_upperstd,
                    h2h_1nearest_LOWERSTD=self.h2h_1nearest_lowerstd,
                    nonh2nonh_1nearest_LOWERSTD=self.nonh2nonh_1nearest_lowerstd,
                    ):
                    return (struct_pmg, "clathrate")
                else:
                    print(f"Generating clathrate failed !!!")
                    return None, None
            else:
                return (struct_pmg, "ramdom structure")
        except:
            print(f"Generating failed !!!")
            return None, None

    # @set_timeout(120, after_timeout)
    def _gen_specify_symbols(self, input_atoms: Atoms):

        spacegroup_number = self.spacegroup_number
        nameofatoms = self.nameofatoms

        try:
            fomula = str(input_atoms.symbols.get_chemical_formula('metal'))
            species_amounts_sites = random.choice(self._group[fomula])
        except KeyError:
            # 为什么捕捉这个错误呢？因为在产生结构的时候很有可能产生一些意料之外的化学配比
            # 例如： 225号空间群使用4a, 4b, 48h 构造晶体结构。
            #   虽然48h 的占位是 [y, y, 0]。如果y=0.5，那么48h就会发生merge，merge成12i
            #   此时225空间群就变成221.
            # 就是因为这样特殊的情况，所以 如果读入的结构的化学配比在 self._group字典 中找到不到的话
            # 就返回一个None
            return None

        amounts = species_amounts_sites[0]
        wyck = species_amounts_sites[1]

        species_radius = list(zip(nameofatoms, self.distancematrix.diagonal()/2.0))

        tm = Tol_matrix()
        for ele_r1, ele_r2 in combinations(species_radius, 2):
            tm.set_tol(ele_r1[0], ele_r2[0], ele_r1[1]+ele_r2[1])

        struc = pyxtal()
        try:
            # print("Now the program will try to create ")
            # print(f"every atoms amouts are {amounts}")
            # print(f"wyckoff positions are {wyck}\n")
            struc.from_random(
                3,
                spacegroup_number,
                nameofatoms,
                amounts,
                factor=2.0,
                sites=wyck,
                tm=tm
            )
            struct_pmg = struc.to_pymatgen()
            if struct_pmg.composition.num_atoms < float(self.maxlimit):
                return struct_pmg
        except Exception as e:
            # print(f"The {fomula} generated by `_gen_specify_symbols` failed ")
            return None

    @classmethod
    def init_from_config(cls, config_d: dict):
        self = cls(
            work_path=config_d["work_path"],
            spacegroup_number=config_d["spacegroup_number"],
            nameofatoms=config_d["nameofatoms"],
            numberofatoms=config_d["numberofatoms"],
            mults_ranges=config_d["mults_ranges"],
            occupied_number_ranges=config_d["occupied_number_ranges"],
            popsize=config_d["popsize"],
            maxlimit=config_d["maxlimit"],
            distancematrix=config_d["distancematrix"],
            clathrate_ratio=config_d["clathrate_ratio"],
            hydrogen_content=config_d['hydrogen_content'],

            remain_H_ratio_upperstd=config_d['remain_H_ratio_upperstd'],
            fraction_of_hydrogen_volume_lowerstd=config_d['fraction_of_hydrogen_volume_lowerstd'],
            shr_num_avg_lowerstd=config_d['shr_num_avg_lowerstd'],
            cage_regularity_avg_upperstd=config_d['cage_regularity_avg_upperstd'],
            h2h_network_regularity_avg_upperstd=config_d['h2h_network_regularity_avg_upperstd'],
            h2h_1nearest_lowerstd=config_d['h2h_1nearest_lowerstd'],
            nonh2nonh_1nearest_lowerstd=config_d['nonh2nonh_1nearest_lowerstd'],
        )
        return self


if __name__ == "__main__":

    # binary
    # binary
    spl_wps = split_wyckoffs(
        work_path='.',
        spacegroup_number = [[229, 229]],
        nameofatoms = ["Ca", "H"],
        numberofatoms = [1, 6],
        mults_ranges=[[1,2],
                      [1,12]],
        occupied_number_ranges=[[1,2],
                               [1,2]],
        popsize=30,
        maxlimit=150,
        distancematrix=[[2.014, 1.908, 1.590],
                        [1.908, 1.802, 1.483],
                        [1.590, 1.483, 1.116],],
        clathrate_ratio=0.0,
        hydrogen_content=0.75
    )
    # ternary
    # spl_wps = split_wyckoffs(
    #     work_path='.',
    #     spacegroup_number = [[195,195]],
    #     nameofatoms = ["Ca", "Y", "La", "H"],
    #     numberofatoms = [1, 1, 1, 6],
    #     mults_ranges=[[1,4],
    #                   [1,4],
    #                   [1,4],
    #                   [4,6]],
    #     occupied_number_ranges=[[1,2],
    #                            [1,2],
    #                            [1,2],
    #                            [1,3]],
    #     popsize=30,
    #     maxlimit=150,
    #     distancematrix=[[2.014, 1.908, 1.590],
    #                     [1.908, 1.802, 1.483],
    #                     [1.590, 1.483, 1.116],],
    #     clathrate_ratio=0.0,
    #     hydrogen_content=0.0,
    # )
    struct, struct_type = spl_wps._gen_randomly()
    print(struct, "\n", struct_type)