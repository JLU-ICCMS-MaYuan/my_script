from argparse import ArgumentParser

from relax import relax
from phono import phono

def set_more_args(parser: ArgumentParser):

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

    # 结构弛豫
    subparsers = parser.add_subparsers(help="subparsers")
    parser_relax = subparsers.add_parser("relax")
    parser_relax.add_argument(
        '-m',
        '--more-argments-about-relax',
        type=str,
        dest='more_args_relax',
        nargs='+',
        help="输入更多关于结构弛豫的参数"
    )
    parser_relax.set_defaults(vasp_workflow=relax)

    # 计算声子谱
    parser_phono = subparsers.add_parser("phono")
    parser_phono.add_argument(
        '-m',
        '--more-argments-about-relax',
        type=str,
        dest='more_args_relax',
        nargs='+',
        help="输入更多关于结构弛豫的参数"
    )
    parser_phono.set_defaults(vasp_workflow=phono) 
    
    args = parser.parse_args()

    return args