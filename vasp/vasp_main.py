#!/public/home/mayuan/miniconda3/envs/cage/bin/python3
#!/work/home/mayuan/miniconda3/envs/cage/bin/python3

import os
import logging
from pathlib import Path
from argparse import ArgumentParser

from vasp_workflow import vasp_workflow


logging.basicConfig(level = logging.INFO,format = '%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

if __name__ == "__main__":
    logger.info("Start qe calculate")
    parser = ArgumentParser(prog="run_vasp")
    # 指明输入文件
    parser.add_argument(
        '-i',
        '--input-file-path',
        type=str,
        default=None,
        dest='input_file_path',
        help="please tell me your POSCAR path, please notice your file format had better to be ***.vasp",
    )
    # 指明压强
    parser.add_argument(
        '-p',
        '--press',
        type=float,
        default=0.0,
        dest='press',
        help="please tell me your press which you were on",
    )
    # 指明工作目录
    parser.add_argument(
        '-w',
        '-work-path',
        type=str,
        default=None,
        dest='work_path',
        help="please tell me your calculated directory, I will put all file there",
    )
    # 指明提交作业的系统
    parser.add_argument(
        '-j',
        '-submit-job-system',
        type=str,
        default="slurm",
        dest='submit_job_system',
        help="please tell me your job submition system, eg: slurm, pbs",
    )

    subparsers_relax = parser.add_subparsers(help="relax struct help")
    parser_relax = subparsers_relax.add_parser("relax", help="relax help")
    parser_relax.add_argument(
        "-e",
        "--encut",
        default=None,
        type=int,
        dest="encut",
        help="请输入截断能"
    )
    parser_relax.add_argument(
        "-k",
        "--kspacing",
        default=None,
        type=float,
        dest="kspacing",
        help="请输入截断能"
    ) 
    parser_relax.add_argument(
        "-rv3",
        "--run-vasp-structure-relex-3steps",
        default=False,
        action="store_true",
        dest="run_vasp_3steps",
        help="是否要运行多步vasp结构优化"
    )
    parser_relax.add_argument(
        "-rvf",
        "--run-vasp-high-fine",
        default=False,
        action="store_true",
        dest="run_vasp_high_fine",
        help="是否要运行一步vasp精细结构优化"
    )
    args = parser.parse_args()
    input_file_path   = args.input_file_path
    press             = args.press
    work_path         = args.work_path
    submit_job_system = args.submit_job_system

    encut              = args.encut
    kspacing           = args.kspacing
    run_vasp_3steps    = args.run_vasp_3steps
    run_vasp_high_fine = args.run_vasp_high_fine

    run_mode_list = [
        "relax",                        # whether run relax.in or not
        "scf",                          # whether run scf.in or not
    ]

    print(
        input_file_path   ,
        press             ,
        work_path         ,
        submit_job_system ,

        encut             ,
        kspacing          ,
        run_vasp_3steps   ,
        run_vasp_high_fine,    
        )

    qe_sc = vasp_workflow(
        input_file_path, 
        work_path, 
        kspacing=kspacing, 
        pressure=press, 
        run_mode=run_mode, 
        submit_job_system=submit_job_system
        )