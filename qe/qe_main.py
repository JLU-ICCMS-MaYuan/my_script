#!/usr/bin/env python

from argparse import ArgumentParser
from argparse import RawTextHelpFormatter

from qe.set_args import set_more_args
from qe.qe_log import qe_log

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

    parser = ArgumentParser(prog="run_qe", formatter_class=RawTextHelpFormatter)
    args = set_more_args(parser)
    qe_log(args.logging_level)
    
    qe = args.qe_workflow(args)
