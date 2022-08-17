#!/public/home/mayuan/miniconda3/envs/cage/bin/python3
#!/work/home/mayuan/miniconda3/envs/cage/bin/python3

'''
vasp_main.py -i ./test/POSCAR -w ./test/ -p 200 -j slurm relax -m encut=400 kspacing=0.3 
vasp_main.py -i ./test/POSCAR -w ./test/ -p 200 -j slurm phono -m encut=400 kspacing=0.3 
'''

import os
import logging
from pathlib import Path
from argparse import ArgumentParser

from set_args import set_more_args

logging.basicConfig(level = logging.INFO,format = '%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

if __name__ == "__main__":

    logger.info("Start vasp calculate")
    
    parser = ArgumentParser(prog="run_vasp")
    args = set_more_args(parser)

    logger.info(f"{args} \n")

    vasp = args.vasp_workflow(args)