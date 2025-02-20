#!/usr/bin/env python
import logging

from argparse import ArgumentParser, RawTextHelpFormatter

from vasp.set_args import set_more_args
from vasp.vasp_logging import vasp_logging
from vasp.vasp_logo import vasp_logo

def pre_words():
    vasp_logo("SCRIPTS 4 VASP")
    print("Start vasp calculate")
    print("This code writed by Ma Yuan, and his e-mail is mayuan@calypso.org.cn\n\n")
    
    
if __name__ == "__main__":

    pre_words()
    
    parser = ArgumentParser(prog="run_vasp", formatter_class=RawTextHelpFormatter)
    args = set_more_args(parser)
    vasp_logging(getattr(logging, args.logging_level.upper()))
    
    vasp = args.vasp_workflow(args)