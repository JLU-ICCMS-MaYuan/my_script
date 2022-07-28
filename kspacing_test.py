#!/public/home/mayuan/miniconda3/envs/cage/bin/python3
'''
测试截断能脚本的灵活性：
    1. 可以灵活的针对不同的 物质体系 进行截断能测
    2. 还没想好
测试截断能脚本的工作思路:
    1. 创建kspacing目录和其中的子目录 0.1 0.2 ... 0.8, 并且创建好INCAR
    2. 将相应的POSCAR POTCAR 放入相应的目录中 kspacing/0.1 kspacing/0.2 ... kspacing/0.8
    3. 进入 kspacing/numbner目录提交任务
Use Example:
    1. 创建kspacing目录, 并准备好POSCAR, POTCAR
        test_kspacing.py -pos ./POSCAR -pot ./POTCAR -krange 0.1 0.9 0.1
    2. 只创建kspacing目录
        test_kspacing.py -krange 0.1 0.9 0.1
    3. 为创建好的kspacing目录准备POSCAR, POTCAR
        kspacing_test.py -pos ./POSCAR -pot ./POTCAR
    4. 完整的计算参数设置
        kspacing_test.py -pos ./500/POSCAR -pot ./500/POTCAR -krange 0.1 0.9 0.1 -rk
        
'''
import os
import shutil
import numpy as np
from argparse import ArgumentParser

def write_incar(incar_dirpath, kspacing):
    incar_filepath = os.path.join(incar_dirpath, "INCAR")
    with open(incar_filepath, "w") as incar:
        incar.write("ISTART   = 0    \n")   
        incar.write("ICHARG   = 2    \n")   
        incar.write("ENCUT    = 600  \n")        
        incar.write("PREC     = A    \n")          
        incar.write("KSPACING = {}   \n".format(str(kspacing))) 
        incar.write("NCORE    = 4    \n")

        incar.write("ISMEAR   = 0    \n")   
        incar.write("SIGMA    = 0.05 \n")   
        incar.write("NELM     = 100  \n")   
        incar.write("NELMIN   = 6    \n")   
        incar.write("EDIFF    = 1e-6 \n")

def slurmkspacing(slurm_dirpath):
    slurm_script_filepath = os.path.join(slurm_dirpath, "slurmKspacing.sh")
    with open(slurm_script_filepath, "w") as slurm:
        slurm.write("#!/bin/sh                           \n")     
        slurm.write("#SBATCH  --job-name=kspacing_test   \n")                         
        slurm.write("#SBATCH  --output=kspacing.out.%j   \n")                       
        slurm.write("#SBATCH  --error=kspacing.err.%j    \n")                      
        slurm.write("#SBATCH  --partition=lhy          \n")     # 注意修改 队列名 lhy lbt               
        slurm.write("#SBATCH  --nodes=1                  \n")             
        slurm.write("#SBATCH  --ntasks=48                \n")               
        slurm.write("#SBATCH  --ntasks-per-node=48       \n")                        
        slurm.write("#SBATCH  --cpus-per-task=1          \n")                     
        slurm.write("\n\n")
        slurm.write("source /work/env/intel2018          \n") 
        slurm.write("ulimit -s unlimited                 \n")
        slurm.write("export I_MPI_ADJUST_REDUCE=3        \n")
        slurm.write("export MPIR_CVAR_COLL_ALIAS_CHECK=0 \n")
        slurm.write("\n\n")
        slurm.write('echo "run fine kspacing test!"\n\n')
        slurm.write('mpirun -n 48 /work/software/vasp.6.1.0/vasp_std  > vasp.log 2>&1 \n')                        


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument(
        '-pos',
        '--poscar-file-path',
        action='store',
        type=str,
        dest='poscar_file',
        help='请输入被测试截断能的结构所在的目录'
    )
    parser.add_argument(
        "-pot",
        "--potcar-file-path",
        default=None,
        type=str,
        dest="potcar_file",
        help="请输入POTCAR的路径"
    )
    parser.add_argument(
        '-krange',
        '--test-kspacing-ranges',
        action='store',
        type=float,
        dest="test_kspacing_ranges",
        nargs='+',
        default=None,
        help="请输入截断能的测试范围, 例如: 0.1 0.9 0.1; "
    )
    parser.add_argument(
        '-rk',
        '--run-vasp-kspacing',
        action='store_true',
        default=False,
        dest="run_vasp_kspacing",
        help="是否进行截断能测试"
    )

    args = parser.parse_args()

    test_kspacing_ranges = args.test_kspacing_ranges
    poscar_file       = args.poscar_file
    potcar_file       = args.potcar_file
    run_vasp_kspacing    = args.run_vasp_kspacing

    # 创建kspacing目录和其中的子目录 0.1 100 ... 1000
    # 并将INCAR写入该目录中
    kspacing_dir = os.path.abspath("kspacing")
    if not os.path.exists(kspacing_dir):
        os.mkdir(kspacing_dir) 
    if test_kspacing_ranges:
        if not test_kspacing_ranges[0] < test_kspacing_ranges[1]:
            test_kspacing_ranges[0], test_kspacing_ranges[1] = test_kspacing_ranges[1], test_kspacing_ranges[0]
        for kp in np.arange(test_kspacing_ranges[0],
                    test_kspacing_ranges[1]+test_kspacing_ranges[2], 
                    test_kspacing_ranges[2]):
            sub_kspacing_dir = os.path.join(kspacing_dir, str(round(kp, ndigits=2)))
            if not os.path.exists(sub_kspacing_dir):
                os.makedirs(sub_kspacing_dir)
            write_incar(sub_kspacing_dir, kp)

    # 将结构POSCAR拷贝进相应的目录
    if poscar_file is not None and potcar_file is not None:
        poscar_file = os.path.join(poscar_file)
        potcar_file = os.path.join(potcar_file)
        for root, dirs, files in os.walk(kspacing_dir):
            if "INCAR" in files:
                sub_kspacing_dir = os.path.abspath(root)
                shutil.copy(poscar_file, sub_kspacing_dir)
                shutil.copy(potcar_file, sub_kspacing_dir)
        
    # 将进入目录提交任务
    if run_vasp_kspacing:
        for root, dirs, files in os.walk(kspacing_dir):
            if  "POSCAR" in files and \
                "POTCAR" in files and \
                "INCAR"  in files:
                sub_kspacing_dir = os.path.abspath(root)
                # 创建提交作业的脚本 并执行
                slurmkspacing(sub_kspacing_dir)
                cwd = os.getcwd()
                os.chdir(root)
                os.system("sbatch slurmKspacing.sh")
                os.chdir(cwd)





