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
from itertools import *
from collections import defaultdict
from pathlib import Path
import signal
from typing import *
from fractions import Fraction

import numpy as np
from pyxtal import pyxtal
from pyxtal.tolerance import Tol_matrix
from pyxtal.symmetry import Group, Wyckoff_position

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


class split_wyckoffs:
    '''
    generator_main.py -w ./Ar-Ne-H-num229-500/ -i ./229.ini method -m mode=specifywps

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
#    nonH_upper_limit = '12e'
    H_lower_limit    = '12d'
    sitesoccupiedrange=[[1,2], 
                        [1,2], 
                        [1,3],] 
    popsize=300 
    maxlimit=150
    distancematrix=[[2.014, 1.908, 1.590],
                    [1.908, 1.802, 1.483],
                    [1.590, 1.483, 1.116],]
    clathrate_ratio=0.75
    '''
    def __init__(
        self              ,
        work_path: Path,
        spacegroup_number: List[List[int]],
        nameofatoms: List[str], 
        numberofatoms: List,
        nonH_upper_mult: int,
        nonH_floor_mult: int,
        H_upper_mult: int,
        H_floor_mult: int,

        sitesoccupiedrange: List[List[int]],
        popsize: int,
        maxlimit: int,
        distancematrix: List[List[float]],
        clathrate_ratio : float,

        hydrogen_content : float
        ):
        
        self.spacegroup_number = spacegroup_number
        self.nameofatoms       = nameofatoms
        self.numberofatoms     = numberofatoms
        self.nonH_upper_mult    = nonH_upper_mult
        self.nonH_floor_mult    = nonH_floor_mult
        self.H_upper_mult       = H_upper_mult
        self.H_floor_mult       = H_floor_mult

        self.sitesoccupiedrange = sitesoccupiedrange

        self.work_path          = Path(work_path)
        if not self.work_path.exists():
            self.work_path.mkdir(parents=True)
        self.popsize            = popsize
        self.maxlimit           = maxlimit
        self.distancematrix     = np.array(distancematrix)
        self.clathrate_ratio    = clathrate_ratio
        self.hydrogen_content   = hydrogen_content

        self.structs = []
        self.satisfied_spgs, self.wyckoffpositions = self.get_spgs(self.spacegroup_number, self.nonH_upper_mult)
        self._group = self.get_group(
            self.satisfied_spgs,
            self.wyckoffpositions,
            self.nameofatoms,
            self.nonH_upper_mult,  
            self.nonH_floor_mult, 
            self.H_upper_mult,     
            self.H_floor_mult,
            self.sitesoccupiedrange,
            )

    def get_spgs(self, spacegroup_number, nonH_upper_mult):
        '''
        判断给定的一系列空间群spacegroup_number的wyckoff position的多重度是否满足输入指定的非氢元素可占据的上限nonH_upper_mult
        satisfied_spgs:  将满足要求的空间群返回
        wyckoffpositions:  将满足要求的空间群的wyckoff position的占据情况做成字典返回
        '''
        allspgs = [spgnum for spgnums in spacegroup_number for spgnum in range(spgnums[0], spgnums[1]+1)]
        wyckoffpositions = defaultdict(dict)
        unsatified_spgs = []
        for spgnum in allspgs:
            spg = Group(spgnum, dim=3)
            wps_list = spg.get_wp_list()

            if nonH_upper_mult >= int(re.search(r'\d+', wps_list[-1]).group()):
                for wpstr in wps_list:
                    wp = spg.get_wyckoff_position(wpstr)
                    # 检查给定的 非氢元素的最高多重度nonH_upper_mult 是否比该空间群下最低多重度还要低
                    if wp.get_dof() > 0:
                        wyckoffpositions[spgnum][wpstr] = True
                    else:
                        wyckoffpositions[spgnum][wpstr] = False
            else:
                unsatified_spgs.append(spgnum)

        print("The multiplicity of wyckoff position of these spacegroup numbers is not satisfy with the requirement that `nonH_upper_mult` < the lowest multiplicity of these spacegroups")           
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
        nonH_upper_mult: int,  
        nonH_floor_mult: int, 
        H_upper_mult: int,     
        H_floor_mult: int,
        sitesoccupiedrange,
        ):
            # check whether the file 'composition.json' exist or not
            import json
            composition_json = self.work_path.joinpath("composition.json")
            read_flag = input("Test the existance of `composition.json`, do you want to read it and not generate it by program(y/Y/yes/Yes)? If you input another words, the program will generate composition.json\n")
            if composition_json.exists() and (read_flag == 'y' or read_flag == 'yes' or read_flag == 'Y' or read_flag == 'Yes'):
                with open(composition_json, 'r') as f:
                    _group = json.load(f)
                    return _group
            else:
                if len(nameofatoms) == 2:
                    _group = self.binary_hydrides(
                        nameofatoms,
                        sitesoccupiedrange,
                        satisfied_spgs,
                        wyckoffpositions,
                        nonH_upper_mult,  
                        nonH_floor_mult, 
                        H_upper_mult,     
                        H_floor_mult,
                    )
                    with open(composition_json, 'w') as f:
                        json.dump(_group, f)
                    return _group
                    # elif len(nameofatoms) == 3:
                    #     _group = self.ternary_hydrides(
                    #         nameofatoms,
                    #         sitesoccupiedrange,
                    #         wyckoffpositions,
                    #         nonH_upper_limit,
                    #         H_lower_limit,
                    #     )
                    #     return _group
                    # elif len(nameofatoms) == 4:
                        # self.quaternary_hydrides(
                        #     nameofatoms,
                        #     sitesoccupiedrange,
                        #     wyckoffpositions,
                        #     nonH_upper_limit,
                        #     H_lower_limit,
                        # )
        
    def binary_hydrides(
        self,
        nameofatoms,
        sitesoccupiedrange,
        satisfied_spgs,
        wyckoffpositions,
        nonH_upper_mult,  
        nonH_floor_mult, 
        H_upper_mult,     
        H_floor_mult,
        ):
        _group = defaultdict(list)
        # 判断 H元素 是否为最后一个元素
        if not nameofatoms[-1] == 'H':
            raise ValueError("you have to set nameofatoms in order as `H element is the last one` !!! Just like nameofatoms=['Ar', 'Ne', 'H']")
        
        nonH1_range = sitesoccupiedrange[:-1][0]
        H_range = sitesoccupiedrange[-1]
        # 考虑该非氢元素可能占据的wps的个数的范围从 `nonH1_range[0]~nonH1_range[-1]+1`
        # 例如： 
        #   nonH1_range = [1, 3]
        # 那么:
        #   range(nonH1_range[0], nonH1_range[-1]+1) = [1, 2, 3]
        for spg in satisfied_spgs:
            wyck_pos = wyckoffpositions[spg]
            for nonH1_num in range(nonH1_range[0], nonH1_range[-1]+1):
                # 考虑非氢元素占据 1, 2, 3个wp占位时，所有可能的排列组合情况
                for nonH1_wps, nonH1_rest in self.get_NonHwps(
                    list(wyck_pos.keys()), 
                    wyck_pos,
                    nonH1_num,
                    nonH_upper_mult,
                    nonH_floor_mult,
                    ):
                    for H_num in range(H_range[0], H_range[-1]+1):
                        for H_wps, H_rest in self.get_Hwps(
                            list(wyck_pos.keys()), 
                            nonH1_rest, 
                            H_num,
                            H_upper_mult,
                            H_floor_mult,
                            ):
                            wyckps  = [nonH1_wps, H_wps]
                            nelems  = self.get_natoms([nonH1_wps, H_wps])
                            # 这里all()函数可以接受一个可迭代对象例如列表元组集合，也可以接受一个生成器表达式，例如左面写的那种
                            # 判断这个可迭代对象中所有元素是否都为 True。如果所有元素都为 True，则返回 True，否则返回 False。   all(i is not None for i in self.numberofatoms)
                            if all(self.numberofatoms) == True:
                                if Fraction(self.numberofatoms[0], self.numberofatoms[1]) == Fraction(nelems[0], nelems[1]):
                                    formula = ''.join(map(str, chain.from_iterable(zip(nameofatoms, nelems))))
                                    formula = Formula(formula).format("metal") # 这里的分子式是没有经过约化的。是多少就是多少，8:8:80 不会变成1:1:10
                                    # 计算氢含量
                                    current_hydrogen_content = nelems[-1] / sum(nelems)
                                    if current_hydrogen_content >= self.hydrogen_content:
                                        _group[formula].append([spg, nelems, wyckps])
                                else:
                                    continue
                            else:
                                formula = ''.join(map(str, chain.from_iterable(zip(nameofatoms, nelems))))
                                formula = Formula(formula).format("metal") # 这里的分子式是没有经过约化的。是多少就是多少，8:8:80 不会变成1:1:10
                                # 计算氢含量
                                current_hydrogen_content = nelems[-1] / sum(nelems)
                                if current_hydrogen_content >= self.hydrogen_content:
                                    _group[formula].append([spg, nelems, wyckps])
        _group = dict(_group)

        if _group:
            return _group
        else:  # 有可能指定的配比在当前空间群和wyckoff组合下并不存在，所以需要提示错误
            msg = f"The specified chemical formula {self.nameofatoms} doesn't exist"
            Formula_Specified_Error(msg)

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
                            # for H_wps, H_rest in self.get_Hwps(nonH2_rest, H_num):
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

                        for nonH3_num in range(nonH3_range[0], nonH3_range[-1]+1):
                            for nonH3_wps, nonH3_rest in self.get_NonHwps(list(wyckoffpositions.keys()), nonH2_rest, nonH3_num, nonH_upper_limit):

                                for H_num in range(H_range[0], H_range[-1]+1):
                                    for H_wps, H_rest in self.get_Hwps(list(wyckoffpositions.keys()), nonH3_rest, H_num, H_lower_limit):
                                    # for H_wps, H_rest in self.get_Hwps(nonH3_rest, H_num):
                                        wyckps  = [nonH1_wps, nonH2_wps, nonH3_wps, H_wps]
                                        nelems  = self.get_natoms([nonH1_wps, nonH2_wps, nonH3_wps, H_wps])
                                        formula = ''.join(map(str, chain.from_iterable(zip(nameofatoms, nelems))))
                                        formula = Formula(formula).format("metal") # 这里的分子式是没有经过约化的。是多少就是多少，8:8:80 不会变成1:1:10
                                        _group[formula].append([nelems, wyckps])
        
        _group = dict(_group)
        return _group

    def get_NonHwps(self, original_wps:List, wps:dict, num:int, nonH_upper_mult:int, nonH_floor_mult:int):
        """
        input parameter:
            original_wps: List. It is consist of keys of `wps`. For example:
                         wyckoff_positions = { '1a': False, '1b' : False, '2c' : True, '4d' : True }
                         original_wps      = ['1a', '1b', '2c', '4d']

            wps: Dict. It is a dictionary of wyckoff positions for non-H specie. For example:
                         wps = {'1a': False, '1b' : False, '2c'} 
                       It is different from `wyckoff_positions`. 
                       The difference is that `wps` is the part of rest of `wyckoff_positions`.
            num: Int.  It determines how many wps will the non-H element eventually occupy 

        yield value:
            yield a combination mode of non-H specie
        """
        allwps = list(wps.keys())
        for comb_list in combinations_with_replacement(allwps, num):
            # print("comb_list", comb_list); input()
            for item in set(comb_list):
                if (wps[item] == False) and (comb_list.count(item) > 1): # False 代表该wps占位只能占据一次
                    # 检查：在将wps排列组合的列表里面, 检查占据情况为False并且存在1次以上
                    # 那么这个wps排列组合的列表就是不正确的
                    break
                if (int(re.search(r'\d+', item).group()) < nonH_floor_mult) or (int(re.search(r'\d+', item).group()) > nonH_upper_mult):
                    # 检查：在将wps排列组合的列表里面, 包含过高 或者or 过低多重度的wps
                    # 那么这个wps排列组合的列表就是不符合要求的
                    break
            else:
                comb_list = sorted(list(comb_list))
                
                # 从输入的字典wps中挑选出被金属原子用于占位的键，并且判定这个键是不是只能被占一次
                # 如果这个键只能被占据一次，那么就将它加入到comb_used这个列表中。表示这些wps只能被金属占据且只能占据一次。
                # 如果这个键能被占据多次(即：value == True), 那么说明它不光可以被金属多次占据，也可以被氢多次占据, 就不把这个wp放入comb_used中。
                comb_used = [key for key, value in wps.items() if (key in comb_list) and (value == False)]
                # 这样全部的wps减去被金属使用过的wps
                restwps = set(allwps) - set(comb_used)
                rest_dict = {key : value for key, value in wps.items() if key in restwps}
                yield comb_list, rest_dict

    def get_Hwps(self, original_wps:List, wps:dict, num:int, H_upper_mult:int, H_floor_mult:int):
        """
        input parameter:
            wps: It is a dictionary of wyckoff positions for non-H specie. For example:
                 wyckoff_positions = { '1a': False, '1b' : False, '2c' : True, '4d' : True }
                 wps = {'2c' : True, '4d' : True }
            num: Int.  It determines how many wps will the H element eventually occupy.
        yield value:
            yield a combination mode of H specie
        
        Please attention:
            When the program permutes the wps of the combination H, 
            it does not take into account `H_lower_limit` of the wps of the H element 
        """
        allwps = list(wps.keys())
        for comb_list in combinations_with_replacement(allwps, num):
            for item in set(comb_list):
                if (wps[item] == False) and (comb_list.count(item) != 1): # False 代表该wps占位只能占据一次
                    break
                if (int(re.search(r'\d+', item).group()) < H_floor_mult) or (int(re.search(r'\d+', item).group()) > H_upper_mult):
                    # 检查：在将wps排列组合的列表里面, 包含过高或者过低多重度的wps
                    # 那么这个wps排列组合的列表就是不符合要求的
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
            return None, None

        # species_radius = list(zip(nameofatoms, self.distancematrix.diagonal()/2.0))
        # tm = Tol_matrix()
        # for ele_r1, ele_r2 in combinations(species_radius, 2):
        #     tm.set_tol(ele_r1[0], ele_r2[0], ele_r1[1]+ele_r2[1])

        struc = pyxtal()
        try:
            logger.info(f"try {spg}-spacegroup {amounts} {wyck}")
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
                if checkclathrate(struct_pmg):
                    return (struct_pmg, "clathrate")
                else:
                    return None, None
            else:
                return (struct_pmg, "ramdom structure")
        except:
            logger.info(f"Generating failed !!!")
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
            struct_pmg = struc.to_pymatgen()
            if struct_pmg.composition.num_atoms < float(self.maxlimit):
                return struct_pmg
        except Exception as e:
            # logger.info(f"The {fomula} generated by `_gen_specify_symbols` failed ")
            return None

    @classmethod
    def init_from_config(cls, config_d: dict):
        self = cls(
            work_path=config_d["work_path"],
            spacegroup_number=config_d["spacegroup_number"],
            nameofatoms=config_d["nameofatoms"],
            numberofatoms=config_d["numberofatoms"],
            nonH_upper_mult=config_d["nonh_upper_mult"], 
            nonH_floor_mult=config_d["nonh_floor_mult"],
            H_upper_mult=config_d["h_upper_mult"],
            H_floor_mult=config_d["h_floor_mult"],
            sitesoccupiedrange=config_d["sitesoccupiedrange"],
            popsize=config_d["popsize"],
            maxlimit=config_d["maxlimit"],
            distancematrix=config_d["distancematrix"],
            clathrate_ratio=config_d["clathrate_ratio"],
            hydrogen_content=config_d['hydrogen_content'],
        )
        return self


if __name__ == "__main__":

    spl_wps = split_wyckoffs(
        work_path='.',
        spacegroup_number = [[2, 230]],
        nameofatoms = ["Ca", "H"],
        numberofatoms = [1, 6],
        nonH_upper_mult = 4,
        nonH_floor_mult = 1,
        H_upper_mult = 80,
        H_floor_mult = 2,
        sitesoccupiedrange=[[1,2],
                            [1,2]],
        popsize=30,
        maxlimit=150,
        distancematrix=[[2.014, 1.908, 1.590],
                        [1.908, 1.802, 1.483],
                        [1.590, 1.483, 1.116],],
        clathrate_ratio=0.0,
        hydrogen_content=0.75
    )
    struct, struct_type = spl_wps._gen_randomly()
    print(struct, struct_type)