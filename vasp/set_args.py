from argparse import ArgumentParser, RawTextHelpFormatter

from vasp.vasp_run import *

def set_more_args(parser: ArgumentParser):

    # 指明输入文件所在的位置
    parser.add_argument(
        '-i',
        '--input-file-path',
        type=str,
        default=None,
        dest='input_file_path',
        help="please tell me your POSCAR path, please notice your file format had better to be ***.vasp\n"
            "attention please:\n"
            "   if you use the `batchrelax` or `batchphono`, your input_file_path have to be a directory.\n"
            "   In the directory, there are your candidate-structures that is waiting to be calculate!!!"
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
            "   1. the directory of `work_path` \n"
            "       when you tell me your work_path, the program will create the directory of `work_path` named:\n"
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
    parser_relax = subparsers.add_parser("relax", formatter_class=RawTextHelpFormatter)
    parser_relax.add_argument(
        '-m',
        '--more-argments-about-relax',
        type=str,
        dest='more_args',
        nargs='+',
        help="输入更多关于结构弛豫的参数, 例如：\n"
            "encut=600 (default)\n"
            "kspacing = 0.3 \n"
            "ismear = 0 \n"
            "sigma = 0.01 \n"
            "ediff = 1e-8 \n"
            "ediffg = -0.01 \n"
            "ibrion = 2 \n"
            "isif = 3 \n"
            "potim = 0.1 \n"
            "nelm = 200 \n"
            "ncore = 1 \n"
            "lreal = Auto \n"
            "mode = None, you can set it to be: rvf rv3 \n"
            "queue = xieyu\n"
    )
    parser_relax.set_defaults(vasp_workflow=vasp_relax)
    
    # 性质计算：自洽，能带，电子态密度
    parser_relax = subparsers.add_parser("eletron", formatter_class=RawTextHelpFormatter)
    parser_relax.add_argument(
        '-m',
        '--more-argments-about-properties',
        type=str,
        dest='more_args',
        nargs='+',
        help="输入更多关于的性质计算：自洽，能带，电子态密度参数, 例如：\n"
            "encut=600 (default)\n"
            "kspacing = 0.3 \n"
            "ismear = 0 \n"
            "sigma = 0.01 \n"
            "ediff = 1e-8 \n"
            "ediffg = -0.01 \n"
            "ibrion = 2 \n"
            "isif = 3 \n"
            "potim = 0.1 \n"
            "nelm = 200 \n"
            "ncore = 1 \n"
            "lreal = Auto \n"
            "mode = None, you can set it to be: rvf rv3 \n"
            "queue = xieyu\n"
    )
    parser_relax.set_defaults(vasp_workflow=vasp_eletron)

    # 计算声子谱
    parser_phono = subparsers.add_parser("phono", formatter_class=RawTextHelpFormatter)
    parser_phono.add_argument(
        '-m',
        '--more-argments-about-relax',
        type=str,
        dest='more_args',
        nargs='+',
        help="输入更多关于声子计算的参数\n"
            "encut=600 (default)\n"
            "kspacing = 0.3  or kpoints=[None, None, None]\n"
            "   The program will prefer to check whether kpoints is set\n"
            "   if there is no kpoints setting, then use kmesh-method to creat kpoints\n"
            "       if you choose to set kspacing, then the program will set kpoints by kspacing! \n"
            "       if you choose to set kpoints, then you need to set kpoints='num1, num2, num3  \n'"
            "supercell='1 1 1\n'"
            "ismear = 0 \n"
            "sigma = 0.01 \n"
            "ediff = 1e-6 when you do phono calculation, you can't to set it!!\n"
            "ediffg = -0.01 when you do phono calculation, you can't to set it!!\n"
            "ibrion = -1 (finit displacement method) or 8 (DFPT method) \n"
            "potim = 0.01 when you do phono calculation, you can't to set it!! It is uesd in displacement-method\n"
            "nelm = 200 \n"
            "lreal = Auto, I suggest you to set it to be 'lreal=.FALSE.' \n"
            "mode = None,  you have `disp` and `dfpt`  \n"
            "queue = xieyu\n"
    )
    parser_phono.set_defaults(vasp_workflow=vasp_phono) 

    # 计算MD
    parser_md = subparsers.add_parser("md", formatter_class=RawTextHelpFormatter)
    parser_md.add_argument(
        '-m',
        '--more-argments-about-md',
        type=str,
        dest='more_args',
        nargs='+',
        help="输入更多关于分子动力学计算的参数\n"
            "POTIM 决定分子动力学的步长, 默认是POTIM=0.5"
    )
    parser_md.set_defaults(vasp_workflow=vasp_phono)

    # 批量结构弛豫
    parser_batch_relax = subparsers.add_parser("batchrelax", formatter_class=RawTextHelpFormatter)
    parser_batch_relax.add_argument(
        '-m',
        '--more-argments-about-batch-relax',
        type=str,
        dest='more_args',
        nargs='+',
        help="输入更多关于批量结构弛豫的参数\n"
            "encut=600 (default)\n"
            "kspacing = 0.3 \n"
            "ismear = 0 \n"
            "sigma = 0.01 \n"
            "ediff = 1e-8 \n"
            "ediffg = -0.01 \n"
            "ibrion = 2 \n"
            "isif = 3 \n"
            "potim = 0.1 \n"
            "nelm = 200 \n"
            "ncore = 1 \n"
            "lreal = Auto \n"
            "mode = None, you can set it to be: rvf rv3 \n"
            "queue = xieyu\n"
    )
    parser_batch_relax.set_defaults(vasp_workflow=vaspbatch_relax) 

    # 批量声子谱
    parser_batch_phono = subparsers.add_parser("batchphono", formatter_class=RawTextHelpFormatter)
    parser_batch_phono.add_argument(
        '-m',
        '--more-argments-about-batch-relax',
        type=str,
        dest='more_args',
        nargs='+',
        help="输入更多关于批量声子计算参数\n"
            "encut=600 (default)\n"
            "kspacing = 0.3  or kpoints=[None, None, None]\n"
            "   The program will prefer to check whether kpoints is set\n"
            "   if there is no kpoints setting, then use kmesh-method to creat kpoints\n"
            "       if you choose to set kspacing, then the program will set kpoints by kspacing! \n"
            "       if you choose to set kpoints, then you need to set kpoints='num1, num2, num3  \n'"
            "supercell='1 1 1\n'"
            "ismear = 0 \n"
            "sigma = 0.01 \n"
            "ediff = 1e-6 when you do phono calculation, you can't to set it!!\n"
            "ediffg = -0.01 when you do phono calculation, you can't to set it!!\n"
            "ibrion = -1 (finit displacement method) or 8 (DFPT method) \n"
            "potim = 0.01 when you do phono calculation, you can't to set it!! It is uesd in displacement-method\n"
            "nelm = 200 \n"
            "lreal = Auto, I suggest you to set it to be 'lreal=.FALSE.' \n"
            "mode = None,  you have `disp` and `dfpt`  \n"
            "queue = xieyu\n"
    )
    parser_batch_phono.set_defaults(vasp_workflow=vaspbatch_phono) 

    # 处理数据 
    parser_data = subparsers.add_parser("data", formatter_class=RawTextHelpFormatter)
    parser_data.add_argument(
        '-m',
        '--more-argments-about-batch-relax',
        type=str,
        dest='more_args',
        nargs='+',
        help="输入更多关于声子计算后处理的的参数\n"
            "supercell has no default, you have to set it by yourself\n"
            "   supercell='num1 num2 num3'\n"
            "   mp='8 8 8'\n"
            "   pdos=AUTO \n"
            "   spectrum=False\n"
            "       It mainly determines whether or not to calculate the phono-spectrum\n"
            "   dos=False\n"
            "       It mainly determines whether or not to calculate the total-dos and project-dos\n"

    )
    parser_data.set_defaults(vasp_workflow=vasp_processdata) 

    # 数据清理
    parser_clear = subparsers.add_parser("clear", formatter_class=RawTextHelpFormatter)
    parser_clear.add_argument(
        '-m',
        '--more-argments-about-clear',
        type=str,
        dest='more_args',
        nargs='+',
        help="输入更多关于清理数据的参数\n"
           "清理全部数据，除了 'POSCAR', 'PPOSCAR', 'POTCAR', 'OUTCAR', 'INCAR*', '*.sh', '*.vasp', '*.slurm'\n"
    )
    parser_clear.set_defaults(vasp_workflow=vasp_clear) 


    args = parser.parse_args()
    return args