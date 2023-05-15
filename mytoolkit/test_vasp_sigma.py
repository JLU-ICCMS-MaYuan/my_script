#!/usr/bin/env python3
import os
import sys
import shutil

import numpy as np

def write_incar(sigma, work_path):
    incar = """# test SIGMA
ISTART   = 0     
ICHARG   = 2     
ENCUT    = 800   
PREC     = A     
SYMPREC  = 1e-05     
ISMEAR   = 0     
SIGMA    = {}    
NELM     = 200   
NELMIN   = 6     
EDIFF    = 1e-8  
EDIFFG   = -0.001   
POTIM    = 0.05
NCORE    = 4     
""".format(sigma)
    
    incar_path = os.path.join(work_path, "INCAR")
    with open(incar_path, "w") as incarfile:
        incarfile.write(incar)

def write_submit(jobtype:str, work_path):

    chaoyin_pbs = """#!/bin/sh
#PBS -N    mayqe                                    
#PBS -q    liuhy         
#PBS -l    nodes=1:ppn=28               
#PBS -j    oe                                      
#PBS -V  
source /public/home/mayuan/intel/oneapi/setvars.sh --force
export I_MPI_ADJUST_REDUCE=3
export MPIR_CVAR_COLL_ALIAS_CHECK=0
ulimit -s unlimited        
cd $PBS_O_WORKDIR    

#killall -9 vasp_std

mpirun -n 28 /public/home/mayuan/software/vasp.6.1.0/bin/vasp_std > vasp.log 2>&1  
"""

    tangB_slurm = """#!/bin/sh
#SBATCH  --job-name=mayqe                      
#SBATCH  --output=log.out                       
#SBATCH  --error=log.err                       
#SBATCH  --partition=lhy          
#SBATCH  --nodes=1                          
#SBATCH  --ntasks=48                          
#SBATCH  --ntasks-per-node=48                          
#SBATCH  --cpus-per-task=1                         

source /work/home/may/intel/oneapi/setvars.sh --force
export I_MPI_ADJUST_REDUCE=3
export MPIR_CVAR_COLL_ALIAS_CHECK=0
ulimit -s unlimited

mpirun -n 48 /work/home/may/software/vasp.6.1.0/bin/vasp_std > vasp.log 2>&1 
"""

    coshare_slurm = """#!/bin/sh
#SBATCH  --job-name=vaspopt-1
#SBATCH  --output=log.out
#SBATCH  --error=log.err
#SBATCH  --partition=normal                       
#SBATCH  --nodes=1                            
#SBATCH  --ntasks=64
#SBATCH  --ntasks-per-node=64
#SBATCH  --cpus-per-task=1                    
#SBATCH  --exclude=node37,node34,node38,node48,node10,node15,node5,node20

source /public/env/mpi_intelmpi-2021.3.0.sh
source /public/env/compiler_intel-compiler-2021.3.0.sh

ulimit -s unlimited
export I_MPI_ADJUST_REDUCE=3
export MPIR_CVAR_COLL_ALIAS_CHECK=0
export I_MPI_FABRICS=shm
export MKL_DEBUG_CPU_TYPE=5

mpirun -n 64 /public/software/apps/vasp/intelmpi/5.4.4/bin/vasp_std > vasp.log 2>&1
"""

    if jobtype == "chaoyin_pbs":
        jobsystem = chaoyin_pbs
    elif jobtype == "tangB_slurm":
        jobsystem = tangB_slurm
    elif jobtype == "coshare_slurm":
        jobsystem = coshare_slurm
    else:
        print("其它机器的任务系统, 你没有准备这个机器的提作业脚本. 程序退出. ")
        sys.exit(1)

    submit_path = os.path.join(work_path, "submit.sh")
    with open(submit_path, "w") as jobfile:
        jobfile.write(jobsystem)

def get_enthalpy_per_atom(outcar_path):
    try:
        dE = float(os.popen(f'grep -s "free  energy   TOTEN" {outcar_path}' + " | tail -n 1 | awk '{print $ 5}'").read().strip('\n'))
        begin_id = os.popen(f'grep -n "position of ions in cartesian coordinates" {outcar_path}').read().split(":")[0]
        N = 0; row_id=int(begin_id)
        while True:
            row_id = row_id+1
            content  = os.popen(f"sed -n '{row_id}p' {outcar_path}").read().strip().split()
            if len(content) == 3:
                N += 1
            else:
                break
        dE_per_atom = dE/N
        return dE_per_atom
    except:
        print("Try to get DeltaH from {} failed So let dE_per_atoms = 1000000000000.0!".format(outcar_path))   
        dE_per_atom = 1000000000000.0
        return dE_per_atom




if __name__ == "__main__":
    
    # jobsystem = "chaoyin_pbs"
    jobsystem = "tangB_slurm"
    # jobsystem = "coshare_slurm"
    print("Note: --------------------")
    print("    你需要在当前目录下准备好: POSCAR, POTCAR")
    print("    测试的SIGMA值分别是: 0.01, 0.02, 0.03, 0.04, 0.05, 0.1, 0.2, 0.5")
    print("    在生成测试脚本时, 默认使用ISMEAR=1  !!!!!!!!! ")
    print("    该脚本不提供自动提任务的命令: 你可以用以下命令提供命令:")
    print("        for i in 0.01 0.02 0.03 0.04 0.05 0.1 0.2 0.5; do cd $i; qsub submit.sh; cd ..; done")
    print("        for i in 0.01 0.02 0.03 0.04 0.05 0.1 0.2 0.5; do cd $i; sbatch submit.sh; cd ..; done")
    print("    为什么要优化 SIGMA 值？")
    print("        若展宽 sigma 太小, 则计算难以收敛；若展宽 sigma 太大, 则会产生多余的熵(entropy), 因此必须选择合适的 sigma 值。")
    print("        ISMEAR 和 SIGMA 这两个关键词要联合起来使用, 前者用来指定smearing 的方法, 后者用来指定 smearing 的展宽")
    print("        ISMEAR 和 SIGMA 的默认值分别为 1 和 0.2。")
    print("        ISMEAR 可能的取值为-5, -4, -3, -2, -1, 0, N (N 表示正整数):")
    print("            ISMEAR=-5, 表示采用 Blochl 修正的四面体方法；")
    print("            ISMEAR=-4, 表示采用四面体方法, 但是没有 Blochl 修正；")
    print("            ISMEAR=-1, 表示采用 Fermi-Dirac smearing 方法")
    print("            ISMEAR= 0, 表示采用 Gaussian smearing 方法")
    print("            ISMEAR= N, 表示采用 Methfessel-Paxton smearing 方法, 其中 N 是表示此方法中的阶数, 一般情况下 N 取 1 或 2, 但是 In most cases and leads to very similarresults 。")
    print("        sigma 值一般在 0.1 - 0.3 eV 范围内。")
    print("        ISMEAR 取值的一些经验：")
    print("            (1) 一般说来, 无论是对何种体系, 进行何种性质的计算, 采用ISMEAR = 0 并选择一个合适的 SIGMA 值, 都能得到合理的结果。")
    print("            (2) 在进行静态计算(能量单点计算, no relaxation in metals ) 或态密度计算且 k 点数目大于 4 时, 取 ISMEAR = -5。")
    print("            (3) 当原胞较大而 k 点数目较小(小于 4 个) 时, 取 ISMEAR = 0, 并选择一个合适的 SIGMA 值。( If the cell is too large (or if you use only a single or two k-points) use ISMEAR=0 in combination with a small SIGMA=0.05)")
    print("            (4) 对半导体或绝缘体, 不论是静态还是结构优化计算, 都取ISMEAR = -5。( Mind: Avoid to use ISMEAR>0 for semiconductors and insulators, since it might cause problems. For insulators use ISMEAR=0 or ISMEAR=-5.)")
    print("            (5) 对金属体系( for relaxations in metals ), 取 ISMEAR = 1 或 2, 并选择一个合适的 SIGMA 值。")
    print("    判断标准")
    print("        熵(entropy)越小越好, 选择 entropy T*S EENTRO 值中最小的那个所对应的 SIGMA 。( SIGMA should be as large as possible keeping the difference between the free energy and the total energy (i.e. the term 'entropy T*S') in the OUTCAR file negligible (1meV/atom). )")
    print("        当 k 点的数目发生变化后, 要重新优化选择 SIGMA 值。")
    print("\nNote: --------------------")
    print("    创建测试VASP的SIGMA输入文件目录以及准备vasp的输入文件")

    sigmas = [0.5, 0.2, 0.1, 0.05, 0.04, 0.03, 0.02, 0.01]
    potcar_path = os.path.abspath("POTCAR")
    poscar_path = os.path.abspath("POSCAR")
    for sigma in sigmas:
        test_path = os.path.abspath(str(sigma))
        if not os.path.exists(test_path):
            os.mkdir(test_path)
        shutil.copy(potcar_path, test_path)
        shutil.copy(poscar_path, test_path)
        write_incar(sigma, test_path)
        write_submit(jobsystem, test_path)

    print("    尝试获得每原子的焓值 dE(meV/atom)")
    sigma_dE = []
    for sigma in sigmas:
        outcar_path = os.path.join(str(sigma), "OUTCAR")
        dE_per_atom = get_enthalpy_per_atom(outcar_path)
        sigma_dE.append([sigma, dE_per_atom])

    sigma_dE = np.array(sigma_dE)
    Ediff = np.diff(sigma_dE, axis=0)[:,-1]
    Ediff = np.insert(Ediff, 0, [0.0])
    Ediff = Ediff[:, np.newaxis]
    sigma_dE_Hdiff = np.hstack((sigma_dE, Ediff))
    if len(sigma_dE_Hdiff) == 8:
        print("{:<12},{:<14},{:<14}".format("sigma", "dE(eV/atom)", "diff(meV/atom)"))
        with open("sigma_dE.csv", 'w') as f:
            f.write("{:<12},{:<14},{:<14}\n".format("sigma", "dE(eV/atom)", "diff(meV/atom)"))
            for sigma, dE_per_atom, Hdiff in sigma_dE_Hdiff:
                f.write("{:<12.3f},{:<14.8f},{:<14.8f}\n".format(sigma, dE_per_atom, Hdiff*1000))
                print("{:<12.3f},{:<14.8f},{:<14.8f}".format(sigma, dE_per_atom, Hdiff*1000))
        print("All OUTCARs are OK, sigma_dE.csv has been wroten in current position")
    else:
        print("If all OUTCARs are OK, sigma_dE.csv will be wroten in current position")

