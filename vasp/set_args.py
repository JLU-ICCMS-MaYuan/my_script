from argparse import ArgumentParser, RawTextHelpFormatter

from vasp_run import vasp_relax, vasp_phono, vaspbatch_relax, vaspbatch_phono, vasp_processdata

def set_more_args(parser: ArgumentParser):

    # 指明输入文件所在的位置
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
        help="please tell me your calculated directory, it will determine two things\n"
            "   1. the directory of `work_underpressure` \n"
            "       when you tell me your work_path, the program will create the directory of `work_underpressure` named:\n"
            "       input_file_name + mode + press\n"
            "           such as: POSCAR-disp-200\n"
            "                    CaH6-dfpt-100\n"
            "                    LaH10-rvf-50\n"
            "   2. the directory of `workpath_pppath` \n"
            "       when you tell me your work_path, the program will create the directory of `workpath_pppath` named:\n"
            "       potcar_lib\n"
            

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
    parser_relax.set_defaults(vasp_workflow=vasp_relax)
    
    # 计算声子谱
    parser_phono = subparsers.add_parser("phono")
    parser_phono.add_argument(
        '-m',
        '--more-argments-about-relax',
        type=str,
        dest='more_args',
        nargs='+',
        help="输入更多关于结构弛豫的参数"
    )
    parser_phono.set_defaults(vasp_workflow=vasp_phono) 

    # 批量结构弛豫
    parser_batch_relax = subparsers.add_parser("batchrelax")
    parser_batch_relax.add_argument(
        '-m',
        '--more-argments-about-batch-relax',
        type=str,
        dest='more_args',
        nargs='+',
        help="输入更多关于批量结构弛豫的参数"
    )
    parser_batch_relax.set_defaults(vasp_workflow=vaspbatch_relax) 

    # 批量声子谱
    parser_batch_relax = subparsers.add_parser("batchphono")
    parser_batch_relax.add_argument(
        '-m',
        '--more-argments-about-batch-relax',
        type=str,
        dest='more_args',
        nargs='+',
        help="输入更多关于批量声子计算参数"
    )
    parser_batch_relax.set_defaults(vasp_workflow=vaspbatch_phono) 

    # 处理数据 
    parser_batch_relax = subparsers.add_parser("data")
    parser_batch_relax.add_argument(
        '-m',
        '--more-argments-about-batch-relax',
        type=str,
        dest='more_args',
        nargs='+',
        help="输入更多关于结构弛豫的参数"
    )
    parser_batch_relax.set_defaults(vasp_workflow=vasp_processdata) 






    args = parser.parse_args()
    return args