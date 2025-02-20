#!/usr/bin/env python
import logging
from argparse import ArgumentParser
from argparse import RawTextHelpFormatter

from qe.set_args import set_more_args
from qe.qe_logging import qe_logging
from qe.qe_logo import qe_logo

def pre_words():
    qe_logo("SCRIPTS 4 QE")
    print("Start qe calculate")
    print("This code writed by Ma Yuan, and his e-mail is mayuan@calypso.org.cn\n\n")


if __name__ == "__main__":
    
    pre_words()

    parser = ArgumentParser(prog="run_qe", formatter_class=RawTextHelpFormatter)
    args = set_more_args(parser)
    qe_logging(getattr(logging, args.logging_level.upper()))
    
    qe = args.qe_workflow(args)
