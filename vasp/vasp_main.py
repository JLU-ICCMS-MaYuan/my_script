#!/usr/bin/env python

'''
vasp_main.py -i ./test/POSCAR          -w ./test/ -j pbs -p 200 relax -m encut=400 kspacing=0.3 mode=rv3
vasp_main.py -i ./test/POSCAR          -w ./test/ -j pbs -p 200 relax -m encut=400 kspacing=0.3 mode=rvf
vasp_main.py -i ./test/POSCAR          -w ./test/ -j pbs -p 200 phono -m supercell=[2,2,2] kpoints=[40,40,40] ismear=1 mode=disp   
vasp_main.py -i ./test/POSCAR          -w ./test/ -j pbs -p 200 phono -m supercell=[2,2,2] kpoints=[40,40,40] ismear=1 mode=dfpt
vasp_main.py -i ./test/0.0/POSCAR-init -w ./test/        -p 200 phono -m mode=dispprog supercell=[2,2,2]
vasp_main.py -i ./test/0.0/POSCAR-init -w ./test/        -p 200 phono -m mode=dfptprog supercell=[2,2,2]
vasp_main.py -i ./test                 -w ./opt   -j pbs -p 200 batchrelax -m encut=400 kspacing=0.3 mode=rv3
vasp_main.py -i ./test                 -w ./opt   -j pbs -p 200 batchrelax -m encut=400 kspacing=0.3 mode=rvf


vasp_main.py -i ./test/POSCAR -w ./test/ -j slurm relax -m encut=400 kspacing=0.3 mode=rv3 queue=lhy
vasp_main.py -i ./test/POSCAR -w ./test/ -j slurm relax -m encut=400 kspacing=0.3 mode=rvf queue=xieyu
vasp_main.py -i ./test/POSCAR -w ./test/ -j slurm phono -m supercell=[2,2,2] kpoints=[40,40,40] mode=disp  
vasp_main.py -i ./test/POSCAR -w ./test/ -j slurm phono -m supercell=[2,2,2] kpoints=[40,40,40] mode=dfpt
vasp_main.py -i ./test/0.0/POSCAR-init -w ./test/ phono -m mode=dispprog supercell=[2,2,2]
vasp_main.py -i ./test/0.0/POSCAR-init -w ./test/ phono -m mode=dfptprog supercell=[2,2,2]

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