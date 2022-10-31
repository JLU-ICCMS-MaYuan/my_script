#!/usr/bin/env python

"""
j=0;x=1; for i in `squeue | awk '{print $1}'`; do  let j+=x; done; echo $j

    qe_main.py -i ./POSCAR          -w ./      -p 200 -j bash relax -m mode=relax-vc core=4 npool=1 queue=local kpoints_dense='24 24 24'
    qe_main.py -i ./200.0/relax.out -w ./200.0 -p 200 -j bash scf   -m mode=scffit   core=4 npool=1 queue=local kpoints_dense='24 24 24'
    qe_main.py -i ./200.0/relax.out -w ./200.0 -p 200 -j bash scf   -m mode=scf      core=4 npool=1 queue=local kpoints_dense='12 12 12'
    # 非自洽计算: 用于计算dos, 使用四面体方法, 可以获得没有展宽的dos结果, 更加准确
    qe_main.py -i ./relax.out -w ./test  -p 压力值 -j pbs scf   -m mode=nscf  core=1 npool=1 queue=local kpoints_dense="20 20 20"

    # no split method
    # 进行不分q点计算
    qe_main.py -i ./200.0/relax.out -w ./200.0 -p 200 -j bash phono -m mode=nosplit core=1 npool=1 queue=local qpoints='4 4 4' dyn0_flag=False 
        # default value : dyn0_flag=False
    # 只计算`*dyn0`
    qe_main.py -i ./200.0/relax.out -w ./200.0 -p 200 -j bash phono -m mode=nosplit core=4 npool=1 queue=local qpoints='4 4 4' dyn0_flag=True
    # spilit method 1
    qe_main.py -i ./200.0/relax.out -w ./200.0/ -p 200 -j bash phono -m mode=split_dyn0 core=1 npool=1 queue=local
    qe_main.py -i ./200.0/relax.out -w ./200.0/ -p 200 -j bash phono -m mode=merge core=1 npool=1 queue=local
    # split method 2
    qe_main.py -i ./200.0/relax.out -w ./200.0/ -p 200 -j bash phono -m mode=split_assignQ core=1 npool=1 queue=local

    # 声子计算其它可使用的参数
    el_ph_nsigma=50  

    qe_main.py -i ./200.0/relax.out -w ./200.0/ -p 200 -j bash  phono -m mode=q2r        core=1 npool=1 queue=local
    qe_main.py -i ./200.0/relax.out -w ./200.0/ -p 200 -j bash  phono -m mode=matdyn     core=1 npool=1 queue=local qinserted=50
    qe_main.py -i ./200.0/relax.out -w ./200.0/ -p 200 -j bash  dos   -m mode=matdyn_dos core=1 npool=1 queue=local qpoints='8 8 8' ndos=500 

    qe_main.py -i ./200.0/relax.out -w ./200.0 -p 200 -j bash   sc    -m mode=McAD       core=1 npool=1 queue=local top_freq=80 deguass=0.5 screen_constant=0.1 smearing_method=1
    qe_main.py -i ./200.0/relax.out -w ./200.0 -p 200 -j bash   sc    -m mode=eliashberg core=1 npool=1 queue=local temperature_points=10000 a2F_dos=a2F.dos3

"""
import logging
from argparse import ArgumentParser
from argparse import RawTextHelpFormatter

from set_args import set_more_args

logging.basicConfig(level = logging.INFO,format = '%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

if __name__ == "__main__":


    logger.info("Start qe calculate")

    parser = ArgumentParser(prog="run_vasp", formatter_class=RawTextHelpFormatter)
    args = set_more_args(parser)

    logger.info(f"{args} \n")
    qe = args.qe_workflow(args)
