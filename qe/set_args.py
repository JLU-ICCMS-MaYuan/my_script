from argparse import ArgumentParser, RawTextHelpFormatter

from qe.qe_run import *

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
    
    subparsers = parser.add_subparsers(help="subparsers")

    # 结构弛豫
    parser_relax = subparsers.add_parser("relax", formatter_class=RawTextHelpFormatter)
    parser_relax.add_argument(
        '-m',
        '--more-argments-about-relax',
        type=str,
        dest='more_args',
        nargs='+',
        help="Input more parameters about structural relaxation, such as:\n"
            "mode=relax-vc (no default) --- The cellular structure relaxes\n"
            "queue=xieyu (default), lhy, lbt --- Submitting a task queue\n"
            "forc_conv_thr=1.0d-5 (default)\n"
            "etot_conv_thr=1.0d-7 (default)\n"
            "smearing=gauss (default)\n"
            "degauss=0.005 (default)\n"
            "ecutwfc=60 (default)\n"
            "ecutrho=720 (default)\n"
            "diagonalization=david (default)\n"
            "conv_thr=1.0d-8 (default)\n"
            "mixing_beta=0.7 (default)\n"
            "press_conv_thr=0.01 (default)\n"
            "kpoints_dense='16 16 16' (default)\n"
            "kpoints_sparse='8 8 8' (default)\n"
    )
    parser_relax.set_defaults(qe_workflow=qe_relax)

    # 自洽迭代
    parser_relax = subparsers.add_parser("scf", formatter_class=RawTextHelpFormatter)
    parser_relax.add_argument(
        '-m',
        '--more-argments-about-relax',
        type=str,
        dest='more_args',
        nargs='+',
        help="Enter more parameters for self-consistent iterations, such as\n"
            "mode=scffit(no default) --- self-consistent calculationself (density kpoints meshing)\n"
            "mode=scf   (no default) --- self-consistent calculationself (sparse kpoints meshing)\n"
            "queue=xieyu(default), lhy, lbt --- Submitting a task queue\n"
            "forc_conv_thr=1.0d-5(default)\n"
            "etot_conv_thr=1.0d-7(default)\n"
            "smearing=gauss(default), but we recommend you to choose 'methfessel-paxton' method\n"
            "degauss=0.005(default)\n"
            "ecutwfc=60(default)\n"
            "ecutrho=720(default)\n"
            "diagonalization=david(default)\n"
            "conv_thr=1.0d-8(default), but we recommend you to set '1.0d-9' when you do scffit or scf\n"
            "mixing_beta=0.7(default), but we recommend you to set '0.8' when you do scffit or scf\n"
            "press_conv_thr=0.01(default)\n"
            "kpoints_dense='16 16 16'(default)\n"
            "kpoints_sparse='8 8 8'(default)\n"
    )
    parser_relax.set_defaults(qe_workflow=qe_scf)

   # 声子计算
    parser_phono = subparsers.add_parser("phono", formatter_class=RawTextHelpFormatter)
    parser_phono.add_argument(
        '-m',
        '--more-argments-about-relax',
        type=str,
        dest='more_args',
        nargs='+',
        help="输入更多关于声子计算的参数\n"
            "tr2_ph = 1.0d-16 (default)\n"
            "el_ph_nsigma=50 (default)\n"
            "el_ph_sigma=0.005 (default)\n"
            "alpha_mix=0.5 (default)\n"
            "dyn0_flag=False (default)"
            "       if you will calculate phono spectrum, you can set dyn0_flag = False. \n"
            "       if you just want to get XXX.dyn0, you can set dyn0_flag = True. \n"
            "qpoints='4 4 4' (default) \n"
            "   attention please, the total number of qpoints had better be 30~40!!!\n"
            "   if you don't know how many q points to set, you can use mytool by command:\n"
            "   tool_main.py -i `x` -w `y` kmesh -m kspacing=0.4 \n"
            "   `x`=`path\POSCAR`\n"
            "   `y`=`the position of directory where KPOINTS will be put`\n"
            "qtot=None (default), you don't need to set this parameter. The program will calculate by itself.\n"
            "qirreduced=None (default), reason is the same just like above\n"
            "qinserted=None (default), reason is the same just like above\n"
            "Now phono part has 4 part !\n"
            "   mode=merge\n"
            "   mode=matdyn, at the moment, you need to set qinserted=50 (default)\n"
            "   mode=nosplit\n"
            "   mode=split_dyn0\n"
            "   mode=split_assignQ\n"
    )
    parser_phono.set_defaults(qe_workflow=qe_phono)

   # eletron计算
    parser_dos = subparsers.add_parser("eletron", formatter_class=RawTextHelpFormatter)
    parser_dos.add_argument(
        '-m',
        '--more-argments-about-eletron',
        type=str,
        dest='more_args',
        nargs='+',
        help="输入电子性质计算的参数\n"
            "kpoints_dense='x x x' (尽量密一点） \n"
            "   `x`=`path\POSCAR`\n"
            "   `y`=`the position of directory where KPOINTS will be put`\n"
            "   mode=eband, eletdos, elepdos\n"
    )
    parser_dos.set_defaults(qe_workflow=qe_eletron)

   # 超导
    parser_phono = subparsers.add_parser("sc",  formatter_class=RawTextHelpFormatter)
    parser_phono.add_argument(
        '-m',
        '--more-argments-about-relax',
        type=str,
        dest='more_args',
        nargs='+',
        help="输入更多关于结构弛豫的参数\n"
            "mode=q2r\n"
            "mode=McAD\n"
            "   if we choose the McMillan Alen-Dynes method to calculate superconducting transition temperature, you need to set: \n"
            "   top_freq=80\n"
            "   deguass=0.5\n"
            "   screen_constant=0.1\n"
            "mode=eliashberg"
            "   if we choose the Eliashberg method to calculate superconducting transition temperature, you need to set: \n"
            "   temperature_steps=10000\n"
            "   a2fdos=True\n"
            "   alpha2fdat=False\n"
    )
    parser_phono.set_defaults(qe_workflow=qe_superconduct)


   # relax scffit scf dyn0file
    parser_prepare = subparsers.add_parser("prepare",  formatter_class=RawTextHelpFormatter)
    parser_prepare.add_argument(
        '-m',
        '--more-argments-about-relax',
        type=str,
        dest='more_args',
        nargs='+',
        help="输入更多关于结构弛豫的参数\n"
            "mode=all"
    )
    parser_prepare.set_defaults(qe_workflow=qe_prepare)




    args = parser.parse_args()
    return args