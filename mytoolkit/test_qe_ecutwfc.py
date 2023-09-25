#!/usr/bin/env python3

import os
import sys
import shutil

import numpy as np

def change_ecutwfc_and_ecutrho(old_scf_in_path, new_scf_in_path, ecutwfc, ecutrho):

    with open(old_scf_in_path, "r") as f:
        lines = f.readlines()

    for idx, line in enumerate(lines):
        if "ecutwfc" in line:
            lines[idx] = f" ecutwfc = {ecutwfc}\n"
        if "ecutrho" in line:
            lines[idx] = f" ecutrho = {ecutrho}\n"
        if "pseudo_dir" in line:
            lines[idx] = " pseudo_dir='../pp'\n"

    with open(new_scf_in_path, "w") as f:
        f.writelines(lines)
    
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
#killall -9 pw.x ph.x

mpirun -np 28 /public/home/mayuan/software/qe-7.1/bin/pw.x -npool 4 <scf.in> scf.out 
"""

    tangB_slurm = """#!/bin/sh 
#SBATCH  --job-name=test                      
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

mpirun -np 48 /work/home/may/software/qe-7.1/bin/pw.x -npool 4 <scf.in> scf.out 
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

mpirun -np 64 /public/software/apps/quantum-espresso/intelmpi/6.7/bin/pw.x -npool 8 <scf.in> scf.out
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

def get_total_energy(scf_out_path):
    """计算并返回每原子的能量 单位的Ry"""
    try:
        Et = float(os.popen(f""" grep -s "!    total energy " {scf_out_path}""" + """  | awk '{print $5}' """).read().strip('\n'))
        N = int(os.popen(f""" grep -s "number of atoms/cell" {scf_out_path}""" + """  | awk '{print $5}' """).read().strip('\n'))
        Et_per_atom = Et/N # 注意单位是Ry
        return Et_per_atom
    except:
        print("Try to get total energy from {} failed So let total energy per atom = 1000000000000.0!".format(scf_out_path))   
        Et_per_atom = 1000000000000.0
        return Et_per_atom


if __name__ == "__main__":
    
    # jobsystem = "chaoyin_pbs"
    jobsystem = "tangB_slurm"
    # jobsystem = "coshare_slurm"

    print("Note: --------------------")
    print("    你需要在当前目录下准备好: scf, pp(pp目录中放好赝势)")
    print("    测试的ecutwfc值分别是: 40  50  60  70  80  90")
    print("    测试的ecutrho值分别是: 480 600 720 840 960 10800")
    print("    该脚本不提供自动提任务的命令: 你可以用以下命令提供命令:")
    print("        for i in 40 50 60 70 80 90; do cd $i; qsub submit.sh;   cd ..; done")
    print("        for i in 40 50 60 70 80 90; do cd $i; sbatch submit.sh; cd ..; done")

    print("Note: --------------------")
    print("    创建测试VASP的ENCUT输入文件目录以及准备vasp的输入文件")
    ecutwfcs = [40, 50, 60, 65, 70, 75, 80, 85, 90, 95, 100, 105, 110]
    old_scf_in_path = os.path.abspath("scf.in")
    for ecutwfc in ecutwfcs:
        test_path = os.path.abspath(str(ecutwfc))
        new_scf_in_path = os.path.join(test_path, "scf.in")
        if not os.path.exists(test_path):
            os.mkdir(test_path)
        change_ecutwfc_and_ecutrho(old_scf_in_path, new_scf_in_path, ecutwfc, ecutwfc*12)
        write_submit(jobsystem, test_path)
    
    print("    尝试获得每原子的焓值 Et(meV/atom)")
    ecutwfc_Et = []
    for ecutwfc in ecutwfcs:
        scf_out_path = os.path.join(str(ecutwfc), "scf.out")
        Et_per_atom = get_total_energy(scf_out_path)
        ecutwfc_Et.append([ecutwfc, Et_per_atom])

    ecutwfc_Et = np.array(ecutwfc_Et)
    # 方式一：通过相邻两个ecutwfc的能量的差判断收敛性
    # Ediff = np.diff(ecutwfc_Et, axis=0)[:,-1]
    # Ediff = np.insert(Ediff, 0, [0.0])
    # Ediff = Ediff[:, np.newaxis]
    # ecutwfc_Et_Ediff = np.hstack((ecutwfc_Et, Ediff))
    
    # 方式二：所有ecutwfc的能量减去最小ecutwcf的能量判断收敛性
    Ediff = ecutwfc_Et[:, 1] - ecutwfc_Et[-1, -1]
    Ediff = Ediff[:, np.newaxis]
    degauss_Et_Ediff = np.hstack((ecutwfc_Et, Ediff))
    ecutwfc_Et_Ediff = np.hstack((ecutwfc_Et, Ediff))
    
    if len(ecutwfc_Et) == 13:
        print("{:<12},{:<14},{:<14},{:<14}".format("ewcutwfc", "Et(Ry/atom)", "Et(eV/atom)", "diff(meV/atom)"))
        with open("ecutwfc.csv", 'w') as f:
            f.write("{:<12},{:<14},{:<14},{:<14}".format("ewcutwfc", "Et(Ry/atom)", "Et(eV/atom)", "diff(meV/atom)"))
            for ewcutwfc, Et_per_atom, Ediff in ecutwfc_Et_Ediff:
                f.write("{:<12.3f},{:<14.8f},{:<14.8f},{:<14.8f}\n".format(ewcutwfc, Et_per_atom, Et_per_atom*13.605693009, Ediff*13.605693009*1000))
                print("{:<12.3f},{:<14.8f},{:<14.8f},{:<14.8f}".format(ewcutwfc, Et_per_atom, Et_per_atom*13.605693009, Ediff*13.605693009*1000))
        print("All scf.outs are OK, ecutwfc.csv has been wroten in current position")
    else:
        print("If all OUTCARs are OK, ecutwfc.csv will be wroten in current position")
