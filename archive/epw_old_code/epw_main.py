#!/usr/bin/env python
import logging
from argparse import ArgumentParser
from argparse import RawTextHelpFormatter

from epw.set_args import set_more_args
from epw.epw_logging import epw_logging
from epw.epw_logo import epw_logo

def pre_words():
    epw_logo("SCRIPTS 4 epw")
    print("Start epw calculate")
    print("This code writed by Ma Yuan, and his e-mail is mayuan@calypso.org.cn\n\n")


if __name__ == "__main__":
    
    pre_words()

    parser = ArgumentParser(prog="run_epw", formatter_class=RawTextHelpFormatter)
    args = set_more_args(parser)
    epw_logging(getattr(logging, args.logging_level.upper()))
    
    epw = args.epw_workflow(args)
