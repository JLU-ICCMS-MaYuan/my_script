#!/usr/bin/env python

import logging
from argparse import ArgumentParser
from argparse import RawTextHelpFormatter

from qe.set_args import set_more_args

logging.basicConfig(level = logging.INFO,format = '%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def pre_words():
    print("--------------------")
    print("--------------------")
    print("--------------------")
    print("--------------------")
    print("Start qe calculate")
    print("This code writed by Ma Yuan, and his e-mail is 1157421359@qq.com")
    print("--------------------")
    print("--------------------")
    print("--------------------")
    print("--------------------")

if __name__ == "__main__":
    
    pre_words()

    parser = ArgumentParser(prog="run_vasp", formatter_class=RawTextHelpFormatter)
    args = set_more_args(parser)

    print(f"{args} \n")
    qe = args.qe_workflow(args)
