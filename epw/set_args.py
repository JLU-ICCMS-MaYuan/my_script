from argparse import ArgumentParser, RawTextHelpFormatter

from epw.epw_run import *

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
        '-w',
        '-work-path',
        type=str,
        default=None,
        dest='work_path',
        help="please tell me your calculated directory\n"
            "   1. if input-file-path is ended with `xxx.vasp`,\n" 
            "       the program will create the directory `xxx/press/`, \n"
            "       the work_path will be work_path/xxx/press/\n"
            "   2. if input-file-path is ended with `relax.out`,\n"
            "       the program will not create any the directory,\n"
            "       the parent path of input_file_path will be the work_path.\n"
            "       such as input-file-path is `home/mayuan/substitute/relax.out`, so the parent path is `home/mayuan/substitute/`\n"
            "       At the moment the work_path still isn't invalid !!! It will determine the position of pp(pseudopotential path) !!!\n"
            "   3. if input-file-path is ended with other formats(such as POSCAR CONTCAR xxx.cif),\n"
            "       the program will only create the directory `/press`, the work_path will be `work_path/press/`\n"
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
    parser.add_argument(
        '-l',
        '-logging-level',
        type=str,
        default="INFO",
        dest='logging_level',
        help="please tell me the output infomation level \n"
             "such as, DEBUG, INFO, WARNING, ERROR, CRITICAL\n"
    ) 
    
    subparsers = parser.add_subparsers(help="subparsers")

    # EPW计算
    parser_epw = subparsers.add_parser("epw_run", formatter_class=RawTextHelpFormatter)
    parser_epw.add_argument(
        '-m',
        '--more-argments-about-epw',
        type=str,
        dest='more_args',
        nargs='+',
        help="输入更多关于epw的参数\n"
    )
    parser_epw.set_defaults(epw_workflow=epw_run)
    
    args = parser.parse_args()
    return args