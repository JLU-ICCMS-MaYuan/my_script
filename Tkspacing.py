#!/usr/bin/env python

import os
import shutil
from argparse import ArgumentParser

import numpy as np

'''
Tkspacing.py -i cellfile/ZnH6.vasp -p POTCAR -j Fopt.sh -r 0.4 0.2 0.05 -run
'''

def write_incar(incar_dirpath, number):
    incar_filepath = os.path.join(incar_dirpath, "INCAR")
    with open(incar_filepath, "w") as incar:
        incar.write("ISTART   = 0    \n")   
        incar.write("ICHARG   = 2    \n")   
        incar.write("ENCUT    = 700  \n")        
        incar.write("PREC     = A    \n")          
        incar.write("KSPACING = {}   \n".format(str(number))) 

        incar.write("ISMEAR   = 0    \n" )   
        incar.write("SIGMA    = 0.01 \n" )   
        incar.write("NELM     = 100  \n" )   
        incar.write("NELMIN   = 6    \n" )   
        incar.write("EDIFF    = 1e-8 \n" )
        incar.write("EDIFFG   = -0.01\n" )
        


if __name__ == "__main__":
    # os.system("encut_test.py -re")
    parser = ArgumentParser()
    parser.add_argument(
        '-i',
        '--directory-of-structure',
        action='store',
        type=str,
        dest='poscar_file',
        help='请输入被测试截断能的结构所在的目录'
    )
    parser.add_argument(
        "-p",
        "--potcar",
        default=None,
        type=str,
        dest="potcar_file",
        help="请输入POTCAR的路径"
    )
    parser.add_argument(
        '-j',
        '--submitjob-path',
        action='store',
        type=str,
        dest="submitjob_path",
        default=None,
        help="提交任务的脚本路径"
    )
    parser.add_argument(
        '-r',
        '--test-ranges',
        action='store',
        type=float,
        dest="test_ranges",
        nargs='+',
        default=None,
        help="请输入测试范围, 例如: 200 1000 50; "
    )
    parser.add_argument(
        '-run',
        '--run-vasp',
        action='store_true',
        default=False,
        dest="run_vasp",
        help="是否进行测试"
    )

    args = parser.parse_args()

    test_ranges   = args.test_ranges
    poscar_file   = args.poscar_file
    potcar_file   = args.potcar_file
    run_vasp      = args.run_vasp
    submitjob_path= args.submitjob_path

    # 创建encut目录和其中的子目录 50 100 ... 1000
    # 并将INCAR写入该目录中
    dst_dir = os.path.abspath("kspacing")
    if not os.path.exists(dst_dir):
        os.mkdir(dst_dir) 
    if test_ranges:
        if not test_ranges[0] < test_ranges[1]:
            test_ranges[0], test_ranges[1] = test_ranges[1], test_ranges[0]
        for number in np.arange(test_ranges[0],
                    test_ranges[1]+test_ranges[2], 
                    test_ranges[2],
                    ):
            sub_dst_dir = os.path.join(dst_dir, str(np.round(number, decimals=4)))
            if not os.path.exists(sub_dst_dir):
                os.makedirs(sub_dst_dir)
            write_incar(sub_dst_dir, number)

    # 将结构POSCAR拷贝进相应的目录, 提交任务命令拷贝到相应目录
    if poscar_file is not None and potcar_file is not None:
        poscar_file = os.path.join(poscar_file)
        potcar_file = os.path.join(potcar_file)
        for root, dirs, files in os.walk(dst_dir):
            if "INCAR" in files:
                sub_dst_dir = os.path.abspath(root)
                os.system(f"cp {poscar_file} {os.path.join(sub_dst_dir, 'POSCAR')}")
                shutil.copy(potcar_file, sub_dst_dir)
                shutil.copy(submitjob_path, sub_dst_dir)

    # 将进入目录提交任务
    if run_vasp:
        for root, dirs, files in os.walk(dst_dir):
            if  "POSCAR" in files and \
                "POTCAR" in files and \
                "INCAR"  in files:
                sub_encut_dir = os.path.abspath(root)
                # 创建提交作业的脚本 并执行
                cwd = os.getcwd()
                os.chdir(root)
                os.system(f"qsub {os.path.split(submitjob_path)[1]}")
                os.chdir(cwd)






