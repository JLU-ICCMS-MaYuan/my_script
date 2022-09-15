#!/usr/bin/env python
import logging
from argparse import ArgumentParser

from set_args import set_more_args

logging.basicConfig(level = logging.INFO,format = '%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

if __name__ == "__main__":

    logger.info("Start toolkit !")
    
    parser = ArgumentParser(prog="toolkit")
    args = set_more_args(parser)

    logger.info(f"{args} \n")

    tool = args.tool_flow(args)
