#!/usr/bin/env python

"""
# check tasks number
j=0;x=1; for i in `squeue | awk '{print $1}'`; do  let j+=x; done; echo $j
pbs
    qe_main.py -i ./CaH6.vasp -w ./test/ -p 压力值 -j pbs relax -m mode=relax-vc 
    qe_main.py -i ./CaH6.vasp -w ./test/ -p 压力值 -j pbs scf   -m mode=scffit
    qe_main.py -i ./CaH6.vasp -w ./test/ -p 压力值 -j pbs scf   -m mode=scf
    # 非自洽计算: 用于计算dos, 使用四面体方法, 可以获得没有展宽的dos结果, 更加准确
    qe_main.py -i ./relax.out -w ./test  -p 压力值 -j pbs scf   -m mode=nscf kpoints_dense=[20,20,20]

    # no split method
    # 进行不分q点计算
    qe_main.py -i ./CaH6.vasp -w ./test/ -p 压力值 -j pbs phono -m mode=nosplit qpoints='4 4 4' dyn0_flag=False 
    # 只计算dyn0_flag
    qe_main.py -w ./CaH6.vasp -w ./test/ -p 压力值 -j pbs phono -m mode=npsplit qpoints='4 4 4' dyn0_flag=True

    # spilit method 1
    qe_main.py -i ./relax.out -w ./ -p 压力值  -j pbs phono -m mode=split_from_dyn0 qpoints='8 8 8'    
    qe_main.py -i ./relax.out -w ./ -p 压力值  -j pbs phono -m mode=merge

    # split method 2
    qe_main.py -i ./relax.out -w ./ -p 压力值  -j pbs phono -m mode=split_specify_q

    # 声子计算其它可使用的参数
    el_ph_nsigma=50  

    qe_main.py -i ./relax.out -w ./ -p 压力值  -j pbs phono -m mode=q2r
    qe_main.py -i ./relax.out -w ./ -p 压力值  -j pbs phono -m mode=matdyn qinserted=50
    qe_main.py -i ./relax.out -w ./ -p 压力值  -j pbs phono -m mode=matdyn_dos

    qe_main.py -i ./relax.out -w ./ -p 压力值  -j pbs superconduct -m mode=McAD       top_freq=80 deguass=0.5 screen_constant=0.1
    qe_main.py -i ./relax.out -w ./ -p 压力值  -j pbs superconduct -m mode=eliashberg temperature_points=10000 a2f_dos=a2F.dos3

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
