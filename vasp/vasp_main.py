#!/public/home/mayuan/miniconda3/envs/cage/bin/python3
#!/work/home/mayuan/miniconda3/envs/cage/bin/python3

'''
qe
'''

import os
import logging
from pathlib import Path
from argparse import ArgumentParser

from vasp_workflow import vasp_workflow
from set_args import set_more_args

logging.basicConfig(level = logging.INFO,format = '%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

if __name__ == "__main__":
    logger.info("Start vasp calculate")
    parser = ArgumentParser(prog="run_vasp")
    args = set_more_args(parser)
    args = parser.parse_args()

    input_file_path   = args.input_file_path
    press             = args.press
    work_path         = args.work_path
    submit_job_system = args.submit_job_system
    more_args_relax   = args.more_args_relax

    run_mode_list = [
        "relax",                        # whether run relax.in or not
        "scf",                          # whether run scf.in or not
    ]

    print(
        input_file_path   ,
        press             ,
        work_path         ,
        submit_job_system ,
        more_args_relax   ,
        vasp_workflow     ,
        )

    #qe_sc = vasp_workflow(
        #input_file_path, 
        #work_path, 
        #pressure=press, 
        #submit_job_system=submit_job_system,
        #more_args_relax=more_args_relax,
        #)