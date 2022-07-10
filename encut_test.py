#!/work/home/mayuan/miniconda3/envs/pyxtal/bin/python3
'''
测试截断能脚本的灵活性：
    1. 可以灵活的针对不同的 物质体系 进行截断能测
    2. 还没想好
测试截断能脚本的工作思路:
    1. 创建encut目录和其中的子目录 50 100 ... 1000, 并且创建好INCAR
    2. 将相应的POSCAR POTCAR 放入相应的目录中 encut/50 encut/100 ... encut/1000
    3. 进入encut/numbner目录提交任务
Use Example:
    1. 创建encut目录, 并准备好POSCAR, POTCAR
        encut_test.py -pos ./POSCAR -pot ./POTCAR -erange 200 1000 50 -re 
    2. 只创建encut目录
        test_encut.py -erange 200 1000 50
    3. 为创建好的encut目录准备POSCAR, POTCAR
        test_encut.py -pos ./POSCAR -pot ./POTCAR
    4. 完整的计算参数设置
        encut_test.py -pos 200-opt/POSCAR -pot 200-opt/POTCAR -erange 200 1000 50 -re
'''
import os
import shutil
from argparse import ArgumentParser


def write_incar(incar_dirpath, encut):
    incar_filepath = os.path.join(incar_dirpath, "INCAR")
    with open(incar_filepath, "w") as incar:
        incar.write("ISTART   = 0    \n")   
        incar.write("ICHARG   = 2    \n")   
        incar.write("ENCUT    = {}   \n".format(str(encut)))        
        incar.write("PREC     = A    \n")          
        incar.write("KSPACING = 0.3  \n") 
        incar.write("NCORE    = 4    \n")

        incar.write("ISMEAR   = 0    \n" )   
        incar.write("SIGMA    = 0.05 \n" )   
        incar.write("NELM     = 100  \n" )   
        incar.write("NELMIN   = 6    \n" )   
        incar.write("EDIFF    = 1e-6 \n" )

def slurmEncut(slurm_dirpath):
    slurm_script_filepath = os.path.join(slurm_dirpath, "slurmEncut.sh")
    with open(slurm_script_filepath, "w") as slurm:
        slurm.write("#!/bin/sh                      \n")     
        slurm.write("#SBATCH  --job-name=encut_test \n")                         
        slurm.write("#SBATCH  --output=encut.out.%j \n")                       
        slurm.write("#SBATCH  --error=encut.err.%j  \n")                      
        slurm.write("#SBATCH  --partition=xieyu     \n")                   
        slurm.write("#SBATCH  --nodes=1             \n")             
        slurm.write("#SBATCH  --ntasks=48           \n")               
        slurm.write("#SBATCH  --ntasks-per-node=48  \n")                        
        slurm.write("#SBATCH  --cpus-per-task=1     \n")                     
        slurm.write("\n\n")
        slurm.write("source /work/env/intel2018\n")
        slurm.write("ulimit -s unlimited                 \n")
        slurm.write("export I_MPI_ADJUST_REDUCE=3        \n")
        slurm.write("export MPIR_CVAR_COLL_ALIAS_CHECK=0 \n")
        slurm.write("\n\n")
        slurm.write('echo "run fine encut test!"\n\n')
        slurm.write('mpirun -n 48 /work/software/vasp.6.1.0/vasp_std  > vasp.log 2>&1 \n')      

if __name__ == "__main__":
    # os.system("encut_test.py -re")
    parser = ArgumentParser()
    parser.add_argument(
        '-pos',
        '--directory-of-structure',
        action='store',
        type=str,
        dest='poscar_file',
        help='请输入被测试截断能的结构所在的目录'
    )
    parser.add_argument(
        "-pot",
        "--potcar",
        default=None,
        type=str,
        dest="potcar_file",
        help="请输入POTCAR的路径"
    )
    parser.add_argument(
        '-erange',
        '--test-encut-ranges',
        action='store',
        type=int,
        dest="test_encut_ranges",
        nargs='+',
        default=None,
        help="请输入截断能的测试范围, 例如: 200 1000 50; "
    )
    parser.add_argument(
        '-re',
        '--run-vasp-encut',
        action='store_true',
        default=False,
        dest="run_vasp_encut",
        help="是否进行截断能测试"
    )

    args = parser.parse_args()

    test_encut_ranges = args.test_encut_ranges
    poscar_file       = args.poscar_file
    potcar_file       = args.potcar_file
    run_vasp_encut    = args.run_vasp_encut

    # 创建encut目录和其中的子目录 50 100 ... 1000
    # 并将INCAR写入该目录中
    encut_dir = os.path.abspath("encut")
    if not os.path.exists(encut_dir):
        os.mkdir(encut_dir) 
    if test_encut_ranges:
        if not test_encut_ranges[0] < test_encut_ranges[1]:
            test_encut_ranges[0], test_encut_ranges[1] = test_encut_ranges[1], test_encut_ranges[0]
        for et in range(test_encut_ranges[0],
                    test_encut_ranges[1]+test_encut_ranges[2], 
                    test_encut_ranges[2]):
            sub_encut_dir = os.path.join(encut_dir, str(et))
            if not os.path.exists(sub_encut_dir):
                os.makedirs(sub_encut_dir)
            write_incar(sub_encut_dir, et)

    # 将结构POSCAR拷贝进相应的目录
    if poscar_file is not None and potcar_file is not None:
        poscar_file = os.path.join(poscar_file)
        potcar_file = os.path.join(potcar_file)
        for root, dirs, files in os.walk(encut_dir):
            if "INCAR" in files:
                sub_encut_dir = os.path.abspath(root)
                shutil.copy(poscar_file, sub_encut_dir)
                shutil.copy(potcar_file, sub_encut_dir)
        
    # 将进入目录提交任务
    if run_vasp_encut:
        for root, dirs, files in os.walk(encut_dir):
            if  "POSCAR" in files and \
                "POTCAR" in files and \
                "INCAR"  in files:
                sub_encut_dir = os.path.abspath(root)
                # 创建提交作业的脚本 并执行
                slurmEncut(sub_encut_dir)
                cwd = os.getcwd()
                os.chdir(root)
                os.system("sbatch slurmEncut.sh")
                os.chdir(cwd)






