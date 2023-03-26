#!/usr/bin/env python

import os
import logging
from pathlib import Path
from argparse import ArgumentParser, RawTextHelpFormatter

from vasp.set_args import set_more_args

logging.basicConfig(level = logging.INFO,format = '%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

if __name__ == "__main__":

    logger.info("Start vasp calculate")
    
    parser = ArgumentParser(prog="run_vasp", formatter_class=RawTextHelpFormatter)
    args = set_more_args(parser)

    vasp = args.vasp_workflow(args)