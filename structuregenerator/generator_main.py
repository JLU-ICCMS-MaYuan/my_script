#!/public/home/mayuan/miniconda3/envs/cage/bin/python3
"""
generator_main.py -w ./Kr-Ne-H-spg229-500/ -i ./input.ini method -m mode=specifywps

input.ini的内容
[base]
spacegroup_number = 229 
nameofatoms = ["Kr", "Ne", "H"] 
optionalsites = [['2a', '6b'], 
                 ['8c','12d'], 
                 ['12e', '16f', '24g', '24h', '48i', '48j', '48k']] 
sitesoccupiedrange=[[1,2], 
                     [1,2], 
                     [1,3],] 
popsize=300 
maxlimit=150
"""
import logging
from argparse import ArgumentParser
from pprint import pprint

from set_args import set_more_args

logging.basicConfig(level = logging.INFO,format = '%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

if __name__ == "__main__":


    logger.info("Start generator structures")

    parser = ArgumentParser(prog="generators")
    args = set_more_args(parser)

    pprint(f"{args} \n")
    generator_res = args.generator(args)