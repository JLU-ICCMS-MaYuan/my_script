from argparse import ArgumentParser

from qe_run import *

def set_more_args(parser: ArgumentParser):

    parser.add_argument(
        '-i',
        '--input-file-path',
        type=str,
        default=None,
        dest='input_file_path',
        help="please tell me your POSCAR path, please notice your file format had better to be ***.vasp.\n"
            "the input-file-path format is free, such as: *.vasp, *cif, POSCAR, CONTCAR, relax.out.... \n"
    )
    parser.add_argument(
        '-p',
        '--press',
        type=float,
        default=0.0,
        dest='press',
        help="please tell me your press where you will calculate",
    )
    parser.add_argument(
        '-w',
        '-work-path',
        type=str,
        default=None,
        dest='work_path',
        help="please tell me your calculated directory\n"
            "   1. if input-file-path is ended with `xxx.vasp`,\n" 
            "       the program will create the directory `xxx/press/`, \n"
            "       the work_underpressure will be work_path/xxx/press/\n"
            "   2. if input-file-path is ended with `relax.out`,\n"
            "       the program will not create any the directory,\n"
            "       the parent path of input_file_path will be the work_underpressure.\n"
            "       such as input-file-path is `home/mayuan/substitute/relax.out`, so the parent path is `home/mayuan/substitute/`\n"
            "       At the moment the work_path is invalid !!! So you don't need to set `-w` parameter\n"
            "   3. if input-file-path is ended with other formats(such as POSCAR CONTCAR xxx.cif),\n"
            "       the program will only create the directory `/press`, the work_underpressure will be `work_path/press/`\n"
    )
    parser.add_argument(
        '-j',
        '-submit-job-system',
        type=str,
        default="slurm",
        dest='submit_job_system',
        help="please tell me your job submition system, \n"
             "such as, slurm, pbs\n"
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

   # 声子计算
    parser_phono = subparsers.add_parser("phono")
    parser_phono.add_argument(
        '-m',
        '--more-argments-about-relax',
        type=str,
        dest='more_args',
        nargs='+',
        help="输入更多关于结构弛豫的参数"
    )
    parser_phono.set_defaults(qe_workflow=qe_phono)

   # 超导
    parser_phono = subparsers.add_parser("superconduct")
    parser_phono.add_argument(
        '-m',
        '--more-argments-about-relax',
        type=str,
        dest='more_args',
        nargs='+',
        help="输入更多关于结构弛豫的参数"
    )
    parser_phono.set_defaults(qe_workflow=qe_superconduct)
    args = parser.parse_args()

    return args