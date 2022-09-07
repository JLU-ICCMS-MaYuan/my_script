#!/usr/bin/env python
# use the python on current python env

# 指定wyckoff position产生结构
"""
generator_main.py -w ./Kr-Ne-H-spg229-500/ -i ./input.ini method -m mode=specifywps

input.ini的内容: 
[specifywps]
spacegroup_number = 229 
nameofatoms = ["Ar", "Ne", "H"] 
optionalsites = [['2a', '6b'], 
                 ['8c','12d'], 
                 ['12e', '16f', '24g', '24h', '48i', '48j', '48k']] 
sitesoccupiedrange=[[1,2], 
                    [1,2], 
                    [1,3],] 
distancematrix=[[],
                [],
                [],]
popsize=300 
maxlimit=150
"""

# 指定初始结构原型进行结构替换
"""
generator_main.py -i ./input.ini -w ./ method -m mode=substitution

input.ini的内容: 
[substitution]
prototype_path = "/public/home/mayuan/code/my_script/test/Mg1B3H20.vasp" # 结构原型的路径
replacement = [["Mg", "Ca", "Sr", "Ba"], # 第一个元素必须是结构原型中的元素
               ["B", "Li", "Na", "K", "Ru", "Cs"]]  # 第一个元素必须是结构原型中的元素
"""

# pso算法演化产生结构
"""
[pso]
nameofatoms = ["Mg", "B", "H"]
popsize=5
distancematrix=[[2.014, 1.908, 1.590],
                [1.908, 1.802, 1.483],
                [1.590, 1.483, 1.116],]

numberOflbest = 4
simthreshold = 0.06
fingerprint = "bcm"
lbest = 1
critic = "enthalpy" 
maxstep= 50
"""

import logging
from argparse import ArgumentParser

from set_args import set_more_args

logging.basicConfig(level = logging.INFO,format = '%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

if __name__ == "__main__":


    logger.info("Start generator structures")

    parser = ArgumentParser(prog="generators")
    args = set_more_args(parser)

    #logger.info(f"{args} \n")
    generator_res = args.generator(args)