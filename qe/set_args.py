from argparse import ArgumentParser

from qe_run import *

def set_more_args(parser: ArgumentParser):

    parser.add_argument(
        '-i',
        '--input-file-path',
        type=str,
        default=None,
        dest='input_file_path',
        help="please tell me your POSCAR path, please notice your file format had better to be ***.vasp",
    )
    parser.add_argument(
        '-p',
        '--press',
        type=float,
        default=0.0,
        dest='press',
        help="please tell me your press which you were on",
    )
    parser.add_argument(
        '-w',
        '-work-path',
        type=str,
        default=None,
        dest='work_path',
        help="please tell me your calculated directory, I will put all file there",
    )
    parser.add_argument(
        '-j',
        '-submit-job-system',
        type=str,
        default="slurm",
        dest='submit_job_system',
        help="please tell me your job submition system, eg: slurm, pbs",
    ) 
    subparsers = parser.add_subparsers(help="subparsers")

    # 结构弛豫
    parser_relax = subparsers.add_parser("relax")
    parser_relax.add_argument(
        '-m',
        '--more-argments-about-relax',
        type=str,
        dest='more_args',
        nargs='+',
        help="输入更多关于结构弛豫的参数"
    )
    parser_relax.set_defaults(qe_workflow=qe_relax)

    # 自洽迭代
    parser_relax = subparsers.add_parser("scf")
    parser_relax.add_argument(
        '-m',
        '--more-argments-about-relax',
        type=str,
        dest='more_args',
        nargs='+',
        help="输入更多关于结构弛豫的参数"
    )
    parser_relax.set_defaults(qe_workflow=qe_scf)
    args = parser.parse_args()

    return args